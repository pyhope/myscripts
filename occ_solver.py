#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import brentq
import argparse

parser = argparse.ArgumentParser(description="Calculate occupancies using Fermi-Dirac distribution")
parser.add_argument("--input", "-i", default="eqp1.dat", help="Input filename")
parser.add_argument("--fermi_level", "-f", type=float, required=True, help="Fermi level (in eV)")
parser.add_argument("--fermi_level_down", "-fd", type=float, default=None, help="Fermi level of spin down component (in eV)")
parser.add_argument("--temperature", "-t", type=float, default=4000.0, help="Temperature")
parser.add_argument("--nelect", "-n", type=int, default=1048, help="Total number of electrons (NELECT in OUTCAR)")
parser.add_argument("-qe", action="store_true", help="Output QE-style eigenvalues and occupancies")
parser.add_argument("-vasp", action="store_true", help="Output VASP-style eigenvalues and occupancies")

args = parser.parse_args()
0
input_filename = args.input
E_F_up = args.fermi_level

if args.fermi_level_down is not None:
    E_F_dn = args.fermi_level_down
else:
    E_F_dn = args.fermi_level

N_total = args.nelect

ev2ry = 1/13.605703976

sigma = args.temperature * 0.000086173303    # Fermi-Dirac smearing width (in eV)

# ---- Fermi-Dirac distribution ----
def fermi_dirac(e, E_F, sigma):
    return 1.0 / (1.0 + np.exp((e - E_F) / sigma))

with open("fermi.dat", "w") as f:
    f.write(f"{E_F_up:.10f}\n{E_F_dn:.10f}\n")

# ---- read input data ----
data = np.loadtxt(input_filename, skiprows=1)

group1 = data[data[:, 0] == 1]
band_index = group1[:, 1].astype(int)
e_up = group1[:, 3]

group2 = data[data[:, 0] == 2]
e_dn = group2[:, 3]

# ---- calculate occupancies ----
occ_up = fermi_dirac(e_up, E_F_up, sigma)
occ_dn = fermi_dirac(e_dn, E_F_dn, sigma)

# ---- calculate total number of electrons in each spin ----
N_up = np.sum(occ_up)
N_dn = np.sum(occ_dn)

# ---- calculate magnetization ----
M = N_up - N_dn

# ---- write output data ----
output_data = np.column_stack((band_index, e_up, e_dn, occ_up, occ_dn))
np.savetxt("eig.dat", output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'], header=f"{N_total} 1 {len(band_index)}\n\n0.0 0.0 0.0 1.0", comments='')

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

# ---- print magnetization ----
print(f"Magnetization M = {M:.10f}")

with open("magnetization.dat", "w") as f:
    f.write(f"{M:.10f}\n")
