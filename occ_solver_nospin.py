#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Calculate occupancies using Fermi-Dirac distribution")
parser.add_argument("--input", "-i", default="eqp1.dat", help="Input filename")
parser.add_argument("--output", "-o", default="eig.dat", help="Output filename")
parser.add_argument("--fermi_level", "-f", type=float, required=True, help="Fermi level (in eV)")
parser.add_argument("--temperature", "-t", type=float, default=4000.0, help="Temperature")
parser.add_argument("-qe", action="store_true", help="Output QE-style eigenvalues and occupancies")
parser.add_argument("-vasp", action="store_true", help="Output VASP-style eigenvalues and occupancies")
parser.add_argument("-orig", action="store_true", help="Use original eigenvalues")

args = parser.parse_args()

input_filename = args.input
E_F = args.fermi_level

ev2ry = 1/13.605703976

sigma = args.temperature * 0.000086173303    # Fermi-Dirac smearing width (in eV)

# ---- Fermi-Dirac distribution ----
def fermi_dirac(e, E_F, sigma):
    return 1.0 / (1.0 + np.exp((e - E_F) / sigma))

# ---- read input data ----
data = np.loadtxt(input_filename, skiprows=1)

group1 = data[data[:, 0] == 1]
band_index = group1[:, 1].astype(int)

if args.orig:
    e = group1[:, 2]
else:
    e = group1[:, 3]

with open("fermi.dat", "w") as f:
    f.write(f"{E_F:.10f}\n")

# ---- calculate occupancies ----
occ = fermi_dirac(e, E_F, sigma)

N = np.sum(occ) * 2

# ---- print magnetization ----
print(f"Total number of electrons = {N:.10f} ({round(N)})")

with open("N.dat", "w") as f:
    f.write(f"{round(N)}\n")

# ---- write output data ----
output_data = np.column_stack((band_index, e, e, occ, occ))
np.savetxt(args.output, output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'], header=f"{round(N)} 1 {len(band_index)}\n\n0.0 0.0 0.0 1.0", comments='')

if args.qe:
    with open("et_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, e):
            f.write(f"{x:<6d}{y*ev2ry:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, e):
            f.write(f"{x:<6d}{y*ev2ry:>20.10f}\n")
    with open("wg_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, occ):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, occ):
            f.write(f"{x:<6d}{y:>20.10f}\n")

if args.vasp:
    with open("eigen_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, e):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, e):
            f.write(f"{x:<6d}{y:>20.10f}\n")
    with open("ferwe_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, occ):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, occ):
            f.write(f"{x:<6d}{y:>20.10f}\n")
