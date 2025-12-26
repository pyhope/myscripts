#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import brentq
import argparse

parser = argparse.ArgumentParser(description="Fermi level solver")
parser.add_argument("--input", "-i", default="eqp1.dat", help="Input filename")
parser.add_argument("--output", "-o", default="eig.dat", help="Output filename")
parser.add_argument("--temperature", "-t", type=float, default=4000.0, help="Temperature")
parser.add_argument("--nelect", "-n", type=float, default=1048.0, help="Total number of electrons (NELECT in OUTCAR)")
parser.add_argument("-qe", action="store_true", help="Output QE-style eigenvalues and occupancies")
parser.add_argument("-vasp", action="store_true", help="Output VASP-style eigenvalues and occupancies")
parser.add_argument("-orig", action="store_true", help="Use original eigenvalues")

args = parser.parse_args()

input_filename = args.input

N_total = args.nelect

ev2ry = 1/13.605703976

sigma = args.temperature * 0.000086173303    # Fermi-Dirac smearing width (in eV)

# ---- Fermi-Dirac distribution ----
def fermi_dirac(e, E_F, sigma):
    return 1.0 / (1.0 + np.exp((e - E_F) / sigma))

# ---- solve for Fermi level ----
def find_fermi_level(eigenvalues, N_target, sigma):
    emin = np.min(eigenvalues) - 5 * sigma
    emax = np.max(eigenvalues) + 5 * sigma
    func = lambda E_F: np.sum(fermi_dirac(eigenvalues, E_F, sigma)) - N_target
    return brentq(func, emin, emax)

# ---- read input data ----
data = np.loadtxt(input_filename, skiprows=1)

group1 = data[data[:, 0] == 1]
band_index = group1[:, 1].astype(int)

if args.orig:
    e = group1[:, 2]
else:
    e = group1[:, 3]

# ---- calculate total number of electrons ----
N = N_total / 2

# ---- find Fermi levels for spin-up and spin-down ----
E_F = find_fermi_level(e, N, sigma)

# ----  calculate occupancies ----
occ = fermi_dirac(e, E_F, sigma)

# ---- write output data ----
output_data = np.column_stack((band_index, e, e, occ, occ))
np.savetxt(args.output, output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'], header=f"{N_total} 1 {len(band_index)}\n\n0.0 0.0 0.0 1.0", comments='')

print(f"Fermi level:   {E_F:.10f} eV")

with open("fermi.dat", "w") as f:
    f.write(f"{E_F:.10f}\n")

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
