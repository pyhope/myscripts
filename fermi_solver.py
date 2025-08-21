#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import brentq
import argparse

parser = argparse.ArgumentParser(description="Fermi level solver")
parser.add_argument("--input", "-i", default="eqp1.dat", help="Input filename")
parser.add_argument("--magnetization", "-m", type=float, required=True, help="Total magnetization")
parser.add_argument("--temperature", "-t", type=float, default=4000.0, help="Temperature")
parser.add_argument("--nelect", "-n", type=float, default=1048.0, help="Total number of electrons (NELECT in OUTCAR)")
parser.add_argument("-qe", action="store_true", help="Output QE-style eigenvalues and occupancies")
parser.add_argument("-vasp", action="store_true", help="Output VASP-style eigenvalues and occupancies")

args = parser.parse_args()

input_filename = args.input

N_total = args.nelect
M = args.magnetization  # total magnetization

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
e_up = group1[:, 3]

group2 = data[data[:, 0] == 2]
e_dn = group2[:, 3]

# ---- calculate total number of electrons ----
N_up = (N_total + M) / 2
N_dn = (N_total - M) / 2

# ---- find Fermi levels for spin-up and spin-down ----
E_F_up = find_fermi_level(e_up, N_up, sigma)
E_F_dn = find_fermi_level(e_dn, N_dn, sigma)

# ----  calculate occupancies ----
occ_up = fermi_dirac(e_up, E_F_up, sigma)
occ_dn = fermi_dirac(e_dn, E_F_dn, sigma)

# ---- write output data ----
output_data = np.column_stack((band_index, e_up, e_dn, occ_up, occ_dn))
np.savetxt("occ.dat", output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'])

print(f"Spin-up Fermi level:   {E_F_up:.10f} eV")
print(f"Spin-down Fermi level: {E_F_dn:.10f} eV")

with open("fermi.dat", "w") as f:
    f.write(f"{E_F_up:.10f}\n{E_F_dn:.10f}\n")

if args.qe:
    with open("et_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, e_up):
            f.write(f"{x:<6d}{y*ev2ry:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, e_dn):
            f.write(f"{x:<6d}{y*ev2ry:>20.10f}\n")
    with open("wg_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, occ_up):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, occ_dn):
            f.write(f"{x:<6d}{y:>20.10f}\n")

if args.vasp:
    with open("eigen_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, e_up):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, e_dn):
            f.write(f"{x:<6d}{y:>20.10f}\n")
    with open("ferwe_in.dat", "w") as f:
        f.write("1\n")
        for x, y in zip(band_index, occ_up):
            f.write(f"{x:<6d}{y:>20.10f}\n")
        f.write("2\n")
        for x, y in zip(band_index, occ_dn):
            f.write(f"{x:<6d}{y:>20.10f}\n")
