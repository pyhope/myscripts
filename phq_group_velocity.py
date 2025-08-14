#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-d", "--dq", type=float, default=0.001, help="q-point spacing")
parser.add_argument("-a", "--a_lat", type=float, default=4.376208478, help="lattice constant a in angstroms")
parser.add_argument("-b", "--b_lat", type=float, default=4.636170450, help="lattice constant b in angstroms")
parser.add_argument("-c", "--c_lat", type=float, default=6.382051746, help="lattice constant c in angstroms")

args = parser.parse_args()

a_lat = args.a_lat * 1e-10
b_lat = args.b_lat * 1e-10
c_lat = args.c_lat * 1e-10

dq = args.dq

def parse_matdyn_modes(filepath):
    q_list = []
    freq_list = []
    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_q = None
    current_freqs = []

    for line in lines:
        if line.strip().startswith("q ="):
            if current_q is not None:
                q_list.append(current_q)
                freq_list.append(current_freqs)
                current_freqs = []
            current_q = tuple(map(float, line.strip().split('=')[1].split()))
        elif "freq" in line and "[THz]" in line:
            match = re.search(r"freq\s+\(\s*\d+\)\s+=\s+([-\d\.]+)", line)
            if match:
                freq_thz = float(match.group(1))
                current_freqs.append(freq_thz)

    if current_q is not None:
        q_list.append(current_q)
        freq_list.append(current_freqs)

    return np.array(q_list), np.array(freq_list)  # (nq, 3), (nq, nmodes)

def compute_velocity(f_plus, f_minus, dq, a_lat):
    df = f_plus - f_minus
    # mask = (f_plus > 0.01) & (f_minus > 0.01)
    v = np.zeros_like(df)
    v = (a_lat * df * 1e12) / (dq * 2)  # Convert to m/s
    return v

tags = ['orig', 'xp', 'xm', 'yp', 'ym', 'zp', 'zm']
freq_data = {}
q_ref = None

for tag in tags:
    path = os.path.join(tag, 'matdyn.modes')
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found")
    q, freqs = parse_matdyn_modes(path)
    freq_data[tag] = freqs
    if tag == 'orig':
        q_ref = q

nq, nmodes = freq_data['orig'].shape

vel_x = compute_velocity(freq_data['xp'], freq_data['xm'], dq, a_lat)
vel_y = compute_velocity(freq_data['yp'], freq_data['ym'], dq, b_lat)
vel_z = compute_velocity(freq_data['zp'], freq_data['zm'], dq, c_lat)

group_vel = np.stack([vel_x, vel_y, vel_z], axis=-1) # m/s

np.save("group_velocity.npy", group_vel)
np.save("frequencies.npy", freq_data['orig'])

with open("group_velocity.txt", "w") as f:
    for i in range(nq):
        f.write(f"# q-point {i+1}: {q_ref[i][0]:.6f} {q_ref[i][1]:.6f} {q_ref[i][2]:.6f}\n")
        f.write("# mode   freq(THz)    vx(m/s)      vy(m/s)      vz(m/s) \n")
        for j in range(nmodes):
            freq0 = freq_data['orig'][i, j] * 33.35641  # Convert THz to cm^-1
            vx, vy, vz = group_vel[i, j]
            f.write(f"{j+1:5d} {freq0:10.4f} {vx:12.2f} {vy:12.2f} {vz:12.2f}\n")
        f.write("\n")
