#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input_filename", type=str, default="matdyn_orig.in",  help="input filename")
parser.add_argument("-d", "--dq", type=float, default=0.001, help="q-point displacement")
args = parser.parse_args()

dq = args.dq
input_file = args.input_filename

with open(input_file, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if line.strip().isdigit():
        nq = int(line.strip())
        q_start = i + 1
        break

q_lines = lines[q_start:q_start + nq]
q_points = np.array([[float(x) for x in line.split()] for line in q_lines])

header = lines[:q_start]
footer = lines[q_start + nq:]

directions = {
    'orig': np.array([0, 0, 0]),
    'xp': np.array([+dq, 0, 0]),
    'xm': np.array([-dq, 0, 0]),
    'yp': np.array([0, +dq, 0]),
    'ym': np.array([0, -dq, 0]),
    'zp': np.array([0, 0, +dq]),
    'zm': np.array([0, 0, -dq]),
}

for label, disp in directions.items():
    new_qs = q_points[:, :3] + disp
    output_file = f"matdyn_{label}.in"
    with open(output_file, "w") as f:
        f.writelines(header)
        for q in new_qs:
            f.write(f"{q[0]:16.9f} {q[1]:16.9f} {q[2]:16.9f}\n")
        f.writelines(footer)
