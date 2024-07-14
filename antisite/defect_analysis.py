#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Analyze the defects and calculate the defect concentration")
parser.add_argument("--input_file", "-i", type=str, default="selected.dump",  help="input .dump file")
parser.add_argument("--atom_num_per_layer", "-n", type=int, default=72,  help="number of atoms per layer in perfect structure")

args = parser.parse_args()
atom_num_per_layer_perfect = args.atom_num_per_layer

def handle_excess_layer(layer, df, dz_dynamic):
    print("Excess layer with %d atoms!" % len(layer))
    sub_layers = []
    current_sub_layer = [layer[0]]
    for idx in layer[1:]:
        if df.loc[idx, "Z"] - df.loc[current_sub_layer[-1], "Z"] > dz_dynamic:
            sub_layers.append(current_sub_layer)
            print("Define a sub-layer with %d atoms" % len(current_sub_layer))
            current_sub_layer = [idx]
        else:
            current_sub_layer.append(idx)
    if current_sub_layer:
        print("Define a sub-layer with %d atoms" % len(current_sub_layer))
        sub_layers.append(current_sub_layer)
    for sub_layer in sub_layers:
        if len(sub_layer) > args.atom_num_per_layer * 1.5:
            return handle_excess_layer(layer, df, dz_dynamic * 0.8)
    return sub_layers

MD_data = mda.Universe(args.input_file, format='LAMMPSDUMP')
MD_data.transfer_to_memory()
frames = MD_data.trajectory
antisite = open("antisite_defect_ratio.txt", "w")
interstitial = open("interstitial_defect_ratio.txt", "w")
info = open("info.txt", "w")
abnormal = open("abnormal_layers.txt", "w")

for f_index, frame in enumerate(frames):
    print("frame %d:" % f_index)

    df = pd.DataFrame({"Index": list(range(1, len(frame) + 1)), "Type": MD_data.atoms.types, "X": frame.positions[:,0], "Y": frame.positions[:,1], "Z": frame.positions[:,2]})
    # df.to_csv("frame_%d.csv" % f_index)
    df = df.loc[(df.Type == "1") | (df.Type == "2")]

    df.sort_values(by=["Z"], inplace=True, ignore_index=True)

    z_np = df.Z.to_numpy()
    dz_np_sorted = np.sort(z_np[1:] - z_np[:-1])
    ddz_np_sorted = dz_np_sorted[1:] - dz_np_sorted[:-1]

    ddz_max_index = ddz_np_sorted.argmax()
    dz_thres = (dz_np_sorted[ddz_max_index + 1] + dz_np_sorted[ddz_max_index] * 1.414) / 2.414
    print("dz_thres:", dz_thres)

    # Handle periodic boundary condition
    z_length = frame.dimensions[2]
    prev_z = None
    for index, row in df.iterrows():
        if index > atom_num_per_layer_perfect and row['Z'] - prev_z > dz_thres:
            df.loc[df['Z'] < row['Z'], 'Z'] += z_length
            break
        prev_z = row['Z']
    df.sort_values(by=["Z"], inplace=True, ignore_index=True)
    
    # Layer definition
    layers = []

    df_iter = df.iterrows()
    current_layer = [next(df_iter)[0]]
    for index, row in df_iter:
        if row.Z - df.loc[current_layer[-1], "Z"] > dz_thres:
            layers.append(current_layer)
            current_layer = [index]
        else:
            current_layer.append(index)
    if current_layer:
        layers.append(current_layer)

    # Handle excess layers
    interstitial_Si_count = 0
    interstitial_Mg_count = 0
    adjusted_layers = []

    for layer in layers:
        atom_num_per_layer = len(layer)
        if atom_num_per_layer < atom_num_per_layer_perfect * 0.4:
            types_in_layer = df.loc[layer, "Type"].values
            interstitial_Si_count += np.sum(types_in_layer == "2")
            interstitial_Mg_count += np.sum(types_in_layer == "1")
            print("Interstitial atoms at Z =" , df.loc[layer[0], "Z"], "with", atom_num_per_layer, "atoms!")
        elif atom_num_per_layer > atom_num_per_layer_perfect * 1.5:
            sub_layers = handle_excess_layer(layer, df, dz_thres * 0.5)
            for sub_layer in sub_layers:
                if len(sub_layer) < atom_num_per_layer_perfect * 0.4:
                    types_in_layer = df.loc[sub_layer, "Type"].values
                    interstitial_Si_count += np.sum(types_in_layer == "2")
                    interstitial_Mg_count += np.sum(types_in_layer == "1")
                    print("Interstitial atoms at Z =" , df.loc[sub_layer[0], "Z"], "with", len(sub_layer), "atoms!")
                else:
                    adjusted_layers.append(sub_layer)
        else:
            adjusted_layers.append(layer)

    interstitial_atoms_count = interstitial_Si_count + interstitial_Mg_count
    print("Interstitial atom count:", interstitial_atoms_count)

    layers = adjusted_layers

    # Calculate antisite defect ratio
    Mg_in_Si_layer = 0
    Si_in_Mg_layer = 0
    total_Si_perfect = 0
    total_Mg_perfect = 0

    print("Total number of layers:", len(layers))
    abnormal_layers = []

    for layer in layers:
        atom_num_per_layer = len(layer)
        if atom_num_per_layer > atom_num_per_layer_perfect * 1.2 or atom_num_per_layer < atom_num_per_layer_perfect * 0.8:
            print("Number of atoms in a layer is abnormal (%d), continue" % (atom_num_per_layer))
            abnormal_df = df.loc[layer]
            abnormal_df.to_csv("frame_%d_atom_%d.csv" % (f_index, atom_num_per_layer))
            abnormal_layers.append(atom_num_per_layer)
            continue
        Si_num = 0
        Mg_num = 0
        for id in layer:
            if df.loc[id, "Type"] == "2":
                Si_num += 1
            elif df.loc[id, "Type"] == "1":
                Mg_num += 1
        if Si_num > Mg_num:
            total_Si_perfect += atom_num_per_layer_perfect
            Mg_in_Si_layer += Mg_num
        else:
            total_Mg_perfect += atom_num_per_layer_perfect
            Si_in_Mg_layer += Si_num

    print("Total number of Si atoms in perfect structure:", total_Si_perfect)
    print("Total number of Mg atoms in perfect structure:", total_Mg_perfect)
    print("# of Mg in Si layer:", Mg_in_Si_layer)
    print("# of Si in Mg layer:", Si_in_Mg_layer)
    print("# of vacancies of Si:", interstitial_Si_count)
    print("# of vacancies of Mg:", interstitial_Mg_count)
    print("Interstitial defect ratio = vacancy defect ratio =", interstitial_atoms_count / (total_Si_perfect + total_Mg_perfect))
    print("Antisite defect ratio:", (Mg_in_Si_layer + Si_in_Mg_layer) / (total_Si_perfect + total_Mg_perfect))
    print()

    info.write("frame %d:\n" % f_index)
    info.write("Total number of layers: %d\n" % len(layers))
    info.write("# of interstitial atoms: %d\n" % interstitial_atoms_count)
    info.write("Abnormal layers: " + str(abnormal_layers) + '\n')
    info.write("# of Si in perfect structure: %d\n" % total_Si_perfect)
    info.write("# of Mg in perfect structure: %d\n" % total_Mg_perfect)
    info.write("# of Mg in Si layer: %d\n" % Mg_in_Si_layer)
    info.write("# of Si in Mg layer: %d\n" % Si_in_Mg_layer)
    info.write("# of vacancies of Si: %d\n" % interstitial_Si_count)
    info.write("# of vacancies of Mg: %d\n" % interstitial_Mg_count)
    info.write("Interstitial defect ratio = vacancy defect ratio = %f\n" % (interstitial_atoms_count / (total_Si_perfect + total_Mg_perfect)))
    info.write("Antisite defect ratio: %f\n" % ((Mg_in_Si_layer + Si_in_Mg_layer) / (total_Si_perfect + total_Mg_perfect)))
    info.write("\n")

    if len(abnormal_layers) > 0:
        abnormal.write('frame %d: ' % f_index + str(abnormal_layers) + '\n')

    antisite.write("%f\n" % ((Mg_in_Si_layer + Si_in_Mg_layer) / (total_Si_perfect + total_Mg_perfect)))
    interstitial.write("%f\n" % (interstitial_atoms_count / (total_Si_perfect + total_Mg_perfect)))

info.close()
antisite.close()
interstitial.close()
abnormal.close()
