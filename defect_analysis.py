#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ase.io.lammpsdata import read_lammps_data
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Calculate the antisite defect concentration")
parser.add_argument("--input_file", "-i", type=str, default="averaged.lmp",  help="input .lmp file")
parser.add_argument("--atom_num_per_layer", "-n", type=int, default=36,  help="number of atoms per layer in perfect structure")
parser.add_argument("--isppv", "-ppv", default=False, action='store_true', help="Defualt: this is brg structure")

args = parser.parse_args()
atom_num_per_layer_perfect = args.atom_num_per_layer
frame = read_lammps_data(args.input_file, style="atomic")

chemical_symbols = []

for atom in frame:
    if atom.symbol == 'H':
        chemical_symbols.append('Mg')
    elif atom.symbol == 'He':
        chemical_symbols.append('Si')
    elif atom.symbol == 'Li':
        chemical_symbols.append('O')

frame.set_chemical_symbols(chemical_symbols)

x_lst = []
y_lst = []
z_lst = []
for atom in frame.positions:
    x_lst.append(atom[0])
    y_lst.append(atom[1])
    z_lst.append(atom[2])

df = pd.DataFrame({"Index": [i for i in range(1, len(frame) + 1)], "Type": frame.get_chemical_symbols(), "X": x_lst, "Y": y_lst, "Z": z_lst})

if args.isppv:
    df = df.loc[((df.Type == "Mg") | (df.Type == "Si"))]
    atom_num_per_layer_perfect = 72
else:
    df = df.loc[(df.X > 12) & ((df.Type == "Mg") | (df.Type == "Si"))]

df.sort_values(by=["Z"], inplace=True, ignore_index=True)
df.to_csv("sorted.csv")

z_np = df.Z.to_numpy()
dz_np_sorted = np.sort(z_np[1:] - z_np[:-1])
ddz_np_sorted = dz_np_sorted[1:] - dz_np_sorted[:-1]

ddz_max_index = ddz_np_sorted.argmax()
dz_thres = (dz_np_sorted[ddz_max_index + 1] + dz_np_sorted[ddz_max_index] * 1.414) / 2.414
print("dz_thres:", dz_thres)

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

# quantities to calculate
si_to_mg = 0
mg_to_si = 0
total_si_perfect = 0
total_mg_perfect = 0

print("Total number of layers:", len(layers))
abnormal_layers = []

with open("layers.txt", "w") as f:
    f.write("Atoms\tSi\tMg\n")
    for layer in layers:
        atom_num_per_layer = len(layer)
        if atom_num_per_layer > atom_num_per_layer_perfect * 1.2 or atom_num_per_layer < atom_num_per_layer_perfect * 0.8:
            print("Number of atoms in a layer is abnormal (%d), continue" % atom_num_per_layer)
            abnormal_layers.append(atom_num_per_layer)
            continue
        si_num = 0
        mg_num = 0
        for id in layer:
            if df.loc[id, "Type"] == "Si":
                si_num += 1
            elif df.loc[id, "Type"] == "Mg":
                mg_num += 1
        if si_num > mg_num:
            total_si_perfect += atom_num_per_layer_perfect
            si_to_mg += mg_num
        else:
            total_mg_perfect += atom_num_per_layer_perfect
            mg_to_si += si_num
        f.write("%d\t%d\t%d\n" % (atom_num_per_layer, si_num, mg_num))
    #print()

print("Total number of Si atoms in perfect struct:", total_si_perfect)
print("Total number of Mg atoms in perfect struct:", total_mg_perfect)
print("Si to Mg:", si_to_mg)
print("Mg to Si:", mg_to_si)
print("antisite defect ratio:", (si_to_mg + mg_to_si) / (total_si_perfect + total_mg_perfect))
with open("info.txt", "w") as f:
    f.write("Total number of layers: %d\n" % len(layers))
    f.write("Abnormal layers: " + str(abnormal_layers) + '\n')
    f.write("Total number of Si atoms in perfect struct: %d\n" % total_si_perfect)
    f.write("Total number of Mg atoms in perfect struct: %d\n" % total_mg_perfect)
    f.write("Si to Mg: %d\n" % si_to_mg)
    f.write("Mg to Si: %d\n" % mg_to_si)
    f.write("antisite defect ratio: %f\n" % ((si_to_mg + mg_to_si) / (total_si_perfect + total_mg_perfect)))

with open("antisite_defect_ratio.txt", "w") as f:
    f.write("%f" % ((si_to_mg + mg_to_si) / (total_si_perfect + total_mg_perfect)))