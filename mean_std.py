#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Compute mean and standard deviation from a specified line to the end of the file.")
parser.add_argument("-i", "--input_file", type=str, default="pressure.dat",  help="input filename")
parser.add_argument("-p", "--plot_file", type=str, default="pressure.png",  help="plot filename")
parser.add_argument("-s", "--start_line", type=int, default="1000", help="Line number to start computing from (1-based indexing).")
args = parser.parse_args()

filename = args.input_file
start_line = args.start_line

with open(filename, 'r') as f:
    lines = f.readlines()[start_line-1:]  # adjust for 0-based indexing

# Convert the lines to a list of floats
values = [float(line.strip()) for line in lines if line.strip()]

# Calculate the mean and standard deviation
mean_value = np.mean(values)
std_value = np.std(values)

plt.subplots(figsize=(8, 8))
values = np.loadtxt(filename)
plt.plot(values)
plt.xlabel('Line Number')
plt.ylabel('Value')
plt.savefig(args.plot_file, dpi=300)

print(f"From line {args.start_line} to the end of the file:")
print(f"Mean: {mean_value:.4f}")
print(f"Standard Deviation: {std_value:.4f}")
