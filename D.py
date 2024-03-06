#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from matplotlib import pyplot as plt
from scipy.stats import linregress

# Argument parser for command line options
parser = argparse.ArgumentParser()
parser.add_argument("--file", "-i", type=str, default='msd_fft.txt', help="input file")
parser.add_argument("--lim1", "-l1", type=float, default=1e5, help="time threshold 1")
parser.add_argument("--lim2", "-l2", type=float, default=5e5, help="time threshold 2")
parser.add_argument("--timestep", "-ts", type=float, default=50, help="time steps")
parser.add_argument("--output", "-o", type=str, default='D_H', help="output filename")
parser.add_argument("--error_output", "-eo", type=str, default='D_err', help="error output filename")

args = parser.parse_args()
timestep = args.timestep

# Function for linear fitting
def fit(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return slope, intercept, r_value**2, std_err

# Function to calculate MSD
def CalcMSD(msd):
    t = np.array(list(range(len(msd)))) * timestep
    index1 = np.argmax(t >= args.lim1)
    index2 = np.argmax(t >= args.lim2)
    
    slope, intercept, r_squared, std_err = fit(t[index1:index2], msd[index1:index2])

    D = slope / 6  # unit: A^2/fs
    D = D * 1e-5   # unit: m^2/s
    D_error = std_err / 6 * 1e-5  # unit: * m^2/s

    print(f"D = %.2e ± %.2e m^2/s" % (D, D_error))
    print(f"R^2 = {r_squared}")
    return D, D_error

# Plotting and analysis
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
msd = np.loadtxt(args.file, skiprows=1)

with open(args.file, "r") as file:
    first_line = file.readline()
    numbers = first_line.lstrip("# ").split()
    number_list = [int(number) for number in numbers]

dataout = open(args.output, "w")
errorout = open(args.error_output, "w")

if msd.ndim == 1:
    ax.plot(np.array(list(range(len(msd)))) * timestep, msd)
    D, D_error = CalcMSD(msd)
    dataout.write(str(D) + ' ')
    errorout.write(str(D_error) + ' ')
elif msd.ndim == 2:
    Data = dict()
    for index, row in enumerate(msd.T):
        ele_num = number_list[index]
        print('Element:', ele_num)
        ax.plot(np.array(list(range(len(row)))) * timestep, row, label='Element: ' + str(ele_num))
        D, D_error = CalcMSD(row)
        Data[ele_num] = [D, D_error]

    for i in range(1, len(Data) + 1):
        dataout.write('%.3e\n' % Data[i][0])
        errorout.write('%.2e\n' % Data[i][1])
    ax.legend()

dataout.close()
errorout.close()

ax.set_xlabel('Time (fs)')
ax.set_ylabel('MSD (A^2)')
plt.savefig('msd.jpg', dpi=300)
