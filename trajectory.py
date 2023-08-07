#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
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
parser.add_argument("--timestep","-ts",type = int, default = 1, help="timestep")

args   = parser.parse_args()
dim = ['x', 'y', 'z']
posx, posy, posz = [], [], []
data = open(args.file, 'r')
line = data.readline()
while line:
    line = data.readline()
    if len(line.split()) >= 3:
            if line.split()[0] == str(args.atom):
                pos_data = line.split()
                posx.append(float(pos_data[2]))
                posy.append(float(pos_data[3]))
                posz.append(float(pos_data[4]))
pos = np.array([posx, posy, posz])


np.savetxt('pos.txt', pos.T, header='    '.join(dim), fmt = '%2.5f')
plt.figure(figsize=(8, 6), dpi=300)

for i in range(len(dim)):
    x = np.array(list(range(len(pos[i]))))*args.timestep
    plt.plot(x, pos[i], label=dim[i])

plt.xlabel('time (fs)')
plt.ylabel('Position (A)')
plt.legend(fancybox=False, edgecolor='black')
plt.savefig('pos.png',bbox_inches='tight')
