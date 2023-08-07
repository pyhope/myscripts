#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Yihang Peng
import numpy as np
import argparse
import glob
#from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str, default = 'msd_fft.txt' ,help="input file")
parser.add_argument("--lim1","-l1",type = float, default = 1e5, help="time threshold 1")
parser.add_argument("--lim2","-l2",type = float, default = 5e5, help="time threshold 2")
parser.add_argument("--timestep","-ts",type = int, default = 50, help="time steps")
parser.add_argument("--output","-of",type = str, default = 'D_H', help="output filename")

args   = parser.parse_args()


def fit(x, y, deg):
    z = np.polyfit(x, y, deg=deg, full=True)
    pn = np.poly1d(z[0])
    R_square = 1 - z[1][0] / sum((y - np.mean(y))** 2)
    return pn, R_square

msd = np.loadtxt(args.file, skiprows=1)
timestep =args.timestep
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
print("D_H = ", D, "* 1e-9 m^2/s")
print("R^2 = ", Rs)

dataout = open(args.output, "w")
dataout.write(str(D)+' ')