#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp1d
import os
import pickle as pkl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--path", "-p", type=str, default='./', help="Path to the directories")
parser.add_argument("--xmin", "-x1", type=float, default=-8, help="Minimum value of x-axis (proximity)")
parser.add_argument("--xmax", "-x2", type=float, default=17, help="Maximum value of x-axis (proximity)")
parser.add_argument("--xnum", "-xn", type=int, default=25001, help="Info[i]ber of data points in x-axis (proximity)")
parser.add_argument("--output_pkl", "-o", type=str, default='data.pkl', help="Output file for the final data in pkl format")

args = parser.parse_args()
path = args.path + '/'

subdirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
int_dirs = []
for subdir in subdirs:
    if subdir.isdigit():
        int_dirs.append(int(subdir))
int_dirs.sort()
if len(int_dirs) == 0:
    print("No directories found")
    exit()

Results = dict()
Info = dict()

for i in int_dirs:
    atomic_num_data = np.loadtxt(path + str(i) + '/sum_proximity_0_0.txt', unpack=True)
    Info[i] = {'O':dict(), 'Mg':dict(), 'Fe':dict(), 'W':dict()}
    ele_list = ['O', 'Mg', 'Fe', 'W']
    if atomic_num_data[7] > atomic_num_data[8]:
        for ele_index, ele in enumerate(ele_list):
            Info[i][ele]['l'], Info[i][ele]['s'], Info[i][ele]['i'] = atomic_num_data[ele_index*3+1], atomic_num_data[ele_index*3+2], atomic_num_data[ele_index*3+3]
    else:
        for ele_index, ele in enumerate(ele_list):
            Info[i][ele]['s'], Info[i][ele]['l'], Info[i][ele]['i'] = atomic_num_data[ele_index*3+1], atomic_num_data[ele_index*3+2], atomic_num_data[ele_index*3+3]
    Info[i]['lw'], Info[i]['chi'] = atomic_num_data[13], atomic_num_data[14]
    if i == int_dirs[0]:
        for ele in ele_list:
            Info[ele] = int(Info[i][ele]['s'] + Info[i][ele]['l'] + Info[i][ele]['i'])

sorted_chi = np.sort([Info[i]['chi'] for i in int_dirs])

int_dirs_selected = []
for i in int_dirs:
    if Info[i]['chi'] > sorted_chi[int(len(int_dirs)*0.2)]:
        int_dirs_selected.append(i)
        atomic_fraction_data = np.loadtxt(path + str(i) + '/atomic_fraction.txt', unpack=True)
        ele_list = ['Mg', 'O', 'Fe', 'W']
        with open(path + str(i) + '/atomic_fraction.txt', 'r') as file:
            lines = file.readlines()
            fields = lines[1].strip().split()
            if len(fields) == 5 and fields[0] == '#':
                ele_list = fields[1:]
            else:
                exit("Error: Invalid format in atomic_fraction.txt")
        
        density_data = np.loadtxt(path + str(i) + '/prox.txt', unpack=True)
        Results[i] = dict()
        proximity, Results[i][ele_list[0]], Results[i][ele_list[1]], Results[i][ele_list[2]], Results[i][ele_list[3]], Results[i]['density'] = atomic_fraction_data[0], atomic_fraction_data[1], atomic_fraction_data[3], atomic_fraction_data[5], atomic_fraction_data[7], density_data[1]
        
        indices = np.where((proximity >= args.xmin) & (proximity <= args.xmax))
        x = proximity[indices]
        f = dict()
        xnew = np.linspace(args.xmin, args.xmax, args.xnum)
        for key in Results[i].keys():
            f[key] = interp1d(x, Results[i][key][indices], kind='linear', fill_value="extrapolate")
            Results[i][key] = f[key](xnew)
    else:
        print(f"Skipping {i} due to low chi value ({Info[i]['chi']})")

average = dict()
for key in Results[int_dirs_selected[0]].keys():
    average[key] = np.zeros(args.xnum)
    for i in int_dirs_selected:
        average[key] += Results[i][key]
    average[key] /= len(int_dirs_selected)

combined_array = np.column_stack((xnew, average['Mg'], average['O'], average['Fe'], average['W'], average['density']))
np.savetxt('average_atomic_fraction.txt', combined_array[::100], fmt='%.4f', header='proximity Mg O Fe W density')

Counts = dict()
for ele in ['Mg', 'O', 'Fe', 'W']:
    Counts[ele] = dict()
    for key in ['s', 'l', 'i']:
        Counts[ele][key] = np.array([Info[i][ele][key] for i in int_dirs_selected])

counts_array = np.column_stack((int_dirs_selected, 
                                Counts['Mg']['s'], Counts['Mg']['l'], Counts['Mg']['i'], 
                                Counts['O']['s'], Counts['O']['l'], Counts['O']['i'], 
                                Counts['Fe']['s'], Counts['Fe']['l'], Counts['Fe']['i'], 
                                Counts['W']['s'], Counts['W']['l'], Counts['W']['i']))
np.savetxt('counts.txt', counts_array, fmt='%d', header='Index Mg_s Mg_l Mg_i O_s O_l O_i Fe_s Fe_l Fe_i W_s W_l W_i')

with open('parameters.txt', 'w') as file:
    file.write('Index lw chi\n')
    for i in int_dirs:
        file.write(f'{i} {Info[i]["lw"]} {Info[i]["chi"]}\n')

with open(path + str(int_dirs[-1]) + '/selected.dump', 'r') as file:
    lines = file.readlines()
    n_atoms = int(lines[3].strip())
    L_x = float(lines[5].split()[1]) - float(lines[5].split()[0])
    L_y = float(lines[6].split()[1]) - float(lines[6].split()[0])
    L_z = float(lines[7].split()[1]) - float(lines[7].split()[0])
    V = L_x * L_y * L_z
    if n_atoms != Info['Mg'] + Info['O'] + Info['Fe'] + Info['W']:
        print("Error: Inconsistent number of atoms")
        exit()

lw = np.mean([Info[i]['lw'] for i in int_dirs_selected])
System = {'n_atoms': n_atoms, 'Mg': Info['Mg'], 'O': Info['O'], 'Fe': Info['Fe'], 'W': Info['W'], 'L_x': L_x, 'L_y': L_y, 'L_z': L_z, 'V': V, 'lw': lw}
with open('system.txt', 'w') as file:
    file.write('n_atoms Mg O Fe W L_x/A L_y/A L_z/A V/A^3 lw/A\n')
    file.write(f'{n_atoms} {Info["Mg"]} {Info["O"]} {Info["Fe"]} {Info["W"]} {L_x:.4f} {L_y:.4f} {L_z:.4f} {V:.1f} {lw:.4f}\n')

final_data = [combined_array.T, Counts, System]
with open(args.output_pkl, 'wb') as file:
    pkl.dump(final_data, file)
