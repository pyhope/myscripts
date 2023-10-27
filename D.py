#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Yihang Peng
import numpy as np
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--file","-i",type=str, default = 'msd_fft.txt' ,help="input file")
parser.add_argument("--lim1","-l1",type = float, default = 1e5, help="time threshold 1")
parser.add_argument("--lim2","-l2",type = float, default = 5e5, help="time threshold 2")
parser.add_argument("--timestep","-ts",type = float, default = 50, help="time steps")
parser.add_argument("--output","-o",type = str, default = 'D_H', help="output filename")

args   = parser.parse_args()
timestep =args.timestep

def fit(x, y, deg):
    z = np.polyfit(x, y, deg=deg, full=True)
    pn = np.poly1d(z[0])
    R_square = 1 - z[1][0] / sum((y - np.mean(y))** 2)
    return pn, R_square

def CalcMSD(msd):
    t = np.array(list(range(len(msd)))) * timestep

    for i in range(len(t)):
        if t[i] >= args.lim1:
            index1 = i
            break
    for i in range(len(t)):
        if t[i] >= args.lim2:
            index2 = i
            break
    index2 = i
    pn, Rs = fit(t[index1:index2], msd[index1:index2], 1)

    D = pn[1]/6 # unit: A^2/fs
    D = 1000* D # unit: A^2/ps
    D = 10 * D # unit: * 1e-9 m^2/s

    print(t[index1:index2])
    print("D = ", D, "* 1e-9 m^2/s")
    print("R^2 = ", Rs)
    return D

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
msd = np.loadtxt(args.file, skiprows=1)
if msd.ndim == 1:
    ax.plot(np.array(list(range(len(msd)))) * timestep, msd)
    D = CalcMSD(msd)
    dataout = open(args.output, "w")
    dataout.write(str(D)+' ')
elif msd.ndim == 2:
    dataout = open(args.output, "w")
    for index, row in enumerate(msd.T):
        print('Element:', index + 1)
        ax.plot(np.array(list(range(len(row)))) * timestep, row, label='Element: ' + str(index + 1))
        D = CalcMSD(row)
        dataout.write(str(D)+'\n')
    ax.legend()

ax.set_xlabel('Time (fs)')
ax.set_ylabel('MSD (A^2)')
plt.savefig('msd.jpg', dpi=300)
