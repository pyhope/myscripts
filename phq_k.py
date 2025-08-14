#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-gv", "--gv_dir", type=str, default="../3_gv", help="Group velocity calculation directory")
parser.add_argument("-t", "--temperature", type=float, default=4000, help="Temperature (K)")
parser.add_argument("-n", "--n_atoms", type=int, default=28, help="Number of atoms in the unit cell")
parser.add_argument("-a", "--a_lat", type=float, default=4.376208478, help="lattice constant a in angstroms")
parser.add_argument("-b", "--b_lat", type=float, default=4.636170450, help="lattice constant b in angstroms")
parser.add_argument("-c", "--c_lat", type=float, default=6.382051746, help="lattice constant c in angstroms")
parser.add_argument("-o", "--output_filename", type=str, default="k.dat", help="output filename")
args = parser.parse_args()

# === Read q-point weights ===
with open(f"{args.gv_dir}/_inputs/q-mesh.dat") as f:
    Nq = int(f.readline().strip())  # First line: number of irreducible q-points
    data = np.loadtxt(f)            # Remaining lines: qx qy qz weight
w_q = data[:, 3]  # q-point weights (already normalized)

# Check normalization
if not np.isclose(np.sum(w_q), 1.0, rtol=1e-6):
    raise ValueError(f"q-point weights are not normalized! Sum = {np.sum(w_q)}")

N_atoms = args.n_atoms  # Number of atoms in the unit cell
Nmodes = 3 * N_atoms  # Number of modes (3N for 3D)
T = args.temperature  # Temperature (K)

# === Lattice constants & volume ===
a_lat = args.a_lat * 1e-10
b_lat = args.b_lat * 1e-10
c_lat = args.c_lat * 1e-10
V = a_lat * b_lat * c_lat  # Volume of the unit cell (m³)

# === Physical constants ===
hbar = 1.054571817e-34  # J·s
kB = 1.380649e-23       # J/K
THz_to_radps = 1e12 * 2 * np.pi  # Convert THz to rad/s

# === Load input data ===
v = np.load(f"{args.gv_dir}/group_velocity.npy")      # (Nq, Nmodes, 3)
freq = np.load(f"{args.gv_dir}/frequencies.npy")      # (Nq, Nmodes), in THz
_, tau = np.loadtxt("tau_fit.dat", unpack=True)
tau = tau.reshape((Nq, Nmodes)) * 1e-12         # Convert to seconds

# === Heat capacity per mode (J/K/m³) ===
omega = freq * THz_to_radps  # rad/s
x = hbar * omega / (kB * T)
x = np.where(x < 1e-6, 1e-6, x)  # Avoid division by zero
C = kB * x**2 * np.exp(x) / (np.exp(x) - 1)**2 / V

# Weighted average total heat capacity
C_total = np.sum(w_q[:, None] * C)
print(f"heat capacity at {T} K: {C_total:.2e} J/K/m³")

# === Sanity check: high-T limit ===
C_highT = 3 * N_atoms * kB / V
print(f"High-T limit heat capacity: {C_highT:.2e} J/K/m³")
print(f"Ratio (C_total / high-T): {C_total / C_highT:.3f}")

# === Thermal conductivity ===
tau = np.where(np.isfinite(tau) & (tau > 0), tau, np.nan)  # mask bad points

assert Nq == len(w_q) == v.shape[0] == freq.shape[0]
assert Nmodes == v.shape[1] == freq.shape[1] == tau.shape[1]

kappa_x = np.nansum(w_q[:, None] * C * v[:, :, 0]**2 * tau)
kappa_y = np.nansum(w_q[:, None] * C * v[:, :, 1]**2 * tau)
kappa_z = np.nansum(w_q[:, None] * C * v[:, :, 2]**2 * tau)
kappa_avg = (kappa_x + kappa_y + kappa_z) / 3
print(f"Thermal conductivity at {T} K:")
print(f"  kappa_x = {kappa_x:.2f} W/m/K")
print(f"  kappa_y = {kappa_y:.2f} W/m/K")
print(f"  kappa_z = {kappa_z:.2f} W/m/K")
print(f"  kappa_avg = {kappa_avg:.2f} W/m/K")

np.savetxt(
    args.output_filename,
    np.array([[kappa_x, kappa_y, kappa_z, kappa_avg]]),
    header="k_x k_y k_z k_avg (W/m/K)",
    fmt="%.4f"
)
