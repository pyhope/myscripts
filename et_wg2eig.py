#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--format", "-f", default="qe", help="Input file format (qe or vasp)")
parser.add_argument("--output", "-o", default="eig.dat", help="Output filename")
parser.add_argument("--total_electrons", "-n", type=float, default=1048, help="Total number of electrons")
args = parser.parse_args()

ry2ev = 13.605703976

def read_spin_file(filename, factor=1.0):
    band_index = []
    values_up = []
    values_dn = []

    with open(filename, "r") as f:
        lines = f.readlines()

    lines = [line.strip() for line in lines if line.strip()]

    mode = None
    for line in lines:
        if line == "1":
            mode = "up"
            continue
        elif line == "2":
            mode = "dn"
            continue

        parts = line.split()
        idx = int(parts[0])
        val = float(parts[1])

        if mode == "up":
            band_index.append(idx)
            values_up.append(val * factor)
        elif mode == "dn":
            values_dn.append(val * factor)

    return band_index, values_up, values_dn 

if args.format == "qe":
    band_index, e_up, e_dn = read_spin_file("et.dat", factor=ry2ev)
    _, occ_up, occ_dn = read_spin_file("wg.dat")
elif args.format == "vasp":
    band_index, e_up, e_dn = read_spin_file("eigen_in.dat")
    _, occ_up, occ_dn = read_spin_file("ferwe_in.dat")
else:
    raise ValueError("Unsupported format. Use 'qe' or 'vasp'.")

print(f"Read {len(band_index)} bands with {len(e_up)} spin-up and {len(e_dn)} spin-down eigenvalues.")

output_data = np.column_stack((band_index, e_up, e_dn, occ_up, occ_dn))
np.savetxt(args.output, output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'], header=f"{round(args.total_electrons)} 1 {len(band_index)}\n\n0.0 0.0 0.0 1.0", comments='')
