#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# author: yhpeng

import numpy as np
import MDAnalysis as mda
import argparse
from matplotlib import rcParams
import matplotlib.pyplot as plt

rcParams['font.size'] = '16'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str, help="input file")
parser.add_argument("--atom","-a",type = int, default = 1, help="The atom to analyze")
parser.add_argument("--format","-ft",type = str, default = 'LAMMPSDUMP', help="file format, e.g., LAMMPSDUMP, PDB")
parser.add_argument("--timestep","-ts",type = int, default = 1, help="timestep")

args   = parser.parse_args()
dim = ['x', 'y', 'z']

def unwrap(idx):
    #Bulow et al., 2020 unwrap method, eqn1
    newcoords = []
    for i in range(len(frames)):
        if i == 0:
            newcoords.append(frames[i].positions[idx])
        else:
            wr_i = frames[i].positions[idx]
            wr_i1 = frames[i-1].positions[idx]
            un_i1 = newcoords[i-1]
            tmp =   np.floor((wr_i - wr_i1)/box +.5)*box
            newcoords.append(un_i1 + wr_i - wr_i1 - tmp)
            
    return np.array(newcoords)

file = args.file
u_md = mda.Universe(file, format=args.format)
u_md.transfer_to_memory()
frames = u_md.trajectory
box    = frames[0].dimensions[0:3]

pos = unwrap(args.atom).T
np.savetxt('pos.txt', pos.T, header='    '.join(dim), fmt = '%2.8f')

plt.figure(figsize=(8, 6), dpi=300)

for i in range(len(dim)):
    x = np.array(list(range(len(pos[i]))))*args.timestep
    plt.plot(x, pos[i], label=dim[i])

plt.xlabel('time (fs)')
plt.ylabel('Position (A)')
plt.legend(fancybox=False, edgecolor='black')
plt.savefig('pos.png',bbox_inches='tight')
