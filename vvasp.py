#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Compute mean and standard deviation from a specified line to the end of the file.")
parser.add_argument("-p", action="store_true", help="Process pressure data")
parser.add_argument("-v", action="store_true", help="Process volume data")
parser.add_argument("-t", action="store_true", help="Process temperature data")
parser.add_argument("-e", action="store_true", help="Process energy data")
parser.add_argument("-s", "--start_line", type=int, default=200, help="Line number to start computing from (1-based indexing).")
args = parser.parse_args()

# Determine the input and output files based on the arguments
if args.p:
    prop = "pressure"
elif args.v:
    prop = "volume"
elif args.t:
    prop = "temperature"
elif args.e:
    prop = "energy"
else:
    raise ValueError("No data type selected. Use -p, -v, or -e.")

start_line = args.start_line
filename = f"{prop}.dat"
plot_filename = f"{prop}.png"

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
plt.xlabel('Time (fs)')
plt.ylabel(prop)
plt.savefig(plot_filename, dpi=300)

print(f"From line {args.start_line} to the end of the file:")
print(f"Mean: {mean_value:.4f}")
if prop == "volume":
    print(f"Mean dimention: {np.cbrt(mean_value):.4f}")
print(f"Standard Deviation: {std_value:.4f}")
