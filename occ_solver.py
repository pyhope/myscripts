#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Calculate occupancies using Fermi-Dirac distribution")
parser.add_argument("--input", "-i", default="eqp1.dat", help="Input filename")
parser.add_argument("--output", "-o", default="eig.dat", help="Output filename")
parser.add_argument("--fermi_level", "-f", type=float, required=True, help="Fermi level (in eV)")
parser.add_argument("--fermi_level_down", "-fd", type=float, default=None, help="Fermi level of spin down component (in eV)")
parser.add_argument("--temperature", "-t", type=float, default=4000.0, help="Temperature")
parser.add_argument("-qe", action="store_true", help="Output QE-style eigenvalues and occupancies")
parser.add_argument("-vasp", action="store_true", help="Output VASP-style eigenvalues and occupancies")
parser.add_argument("-orig", action="store_true", help="Use original eigenvalues")
parser.add_argument("--midgap", "-mg", action="store_true",
                    help="Re-center E_F_up/E_F_dn to the mean of VBM/CBM found in e_up/e_dn relative to current E_Fs")

args = parser.parse_args()

input_filename = args.input
E_F_up = args.fermi_level

E_F_dn = args.fermi_level_down if args.fermi_level_down is not None else args.fermi_level

ev2ry = 1/13.605703976

sigma = args.temperature * 0.000086173303    # Fermi-Dirac smearing width (in eV)

# ---- Fermi-Dirac distribution ----
def fermi_dirac(e, E_F, sigma):
    return 1.0 / (1.0 + np.exp((e - E_F) / sigma))

def find_vbm_cbm_relative_to_ef(evals: np.ndarray, ef: float):
    below = evals[evals <= ef]
    above = evals[evals >= ef]
    vbm = below.max() if below.size > 0 else evals.min()
    cbm = above.min() if above.size > 0 else evals.max()
    if cbm < vbm:
        vbm, cbm = cbm, vbm
    return vbm, cbm

# ---- read input data ----
data = np.loadtxt(input_filename, skiprows=1)

group1 = data[data[:, 0] == 1]
band_index = group1[:, 1].astype(int)
group2 = data[data[:, 0] == 2]

if args.orig:
    e_up = group1[:, 2]
    e_dn = group2[:, 2]
else:
    e_up = group1[:, 3]
    e_dn = group2[:, 3]

if args.midgap:
    vbm_up, cbm_up = find_vbm_cbm_relative_to_ef(e_up, E_F_up)
    new_E_F_up = 0.5 * (vbm_up + cbm_up)

    vbm_dn, cbm_dn = find_vbm_cbm_relative_to_ef(e_dn, E_F_dn)
    new_E_F_dn = 0.5 * (vbm_dn + cbm_dn)

    print("=== Midgap Fermi-level recentering enabled ===")
    print(f"[UP] E_F (old) = {E_F_up:.10f}, E_F (new) = {new_E_F_up:.10f}")
    print(f"[DN] E_F (old) = {E_F_dn:.10f}, E_F (new) = {new_E_F_dn:.10f}")

    E_F_up = new_E_F_up
    E_F_dn = new_E_F_dn

with open("fermi.dat", "w") as f:
    f.write(f"{E_F_up:.10f}\n{E_F_dn:.10f}\n")

# ---- calculate occupancies ----
occ_up = fermi_dirac(e_up, E_F_up, sigma)
occ_dn = fermi_dirac(e_dn, E_F_dn, sigma)

# ---- calculate total number of electrons in each spin ----
N_up = np.sum(occ_up)
N_dn = np.sum(occ_dn)

# ---- calculate magnetization ----
M = N_up - N_dn
N_total = N_up + N_dn

# ---- print magnetization ----
print(f"Magnetization M = {M:.10f}")
print(f"Total number of electrons = {N_total:.10f} ({round(N_total)})")

with open("M_N.dat", "w") as f:
    f.write(f"{M:.10f}\n")
    f.write(f"{round(N_total)}\n")

# ---- write output data ----
output_data = np.column_stack((band_index, e_up, e_dn, occ_up, occ_dn))
np.savetxt(args.output, output_data, fmt=['%d', '%.10f', '%.10f', '%.10f', '%.10f'], header=f"{round(N_total)} 1 {len(band_index)}\n\n0.0 0.0 0.0 1.0", comments='')

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
