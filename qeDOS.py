#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import glob
import os
import shutil

# === Create output directory ===
output_dir = "dos_total"
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory created: {output_dir}")

# === Read Fermi energy ===
with open("Ef.dat", "r") as f:
    fermi_energy = float(f.readline().strip())
print(f"Fermi energy read: {fermi_energy:.6f} eV")

# === Copy Ef.dat to output directory ===
shutil.copy("Ef.dat", os.path.join(output_dir, "Ef.dat"))

# === Initialize containers ===
element_orbital_sum = {}  # e.g. Fe_d
element_total_sum = {}    # e.g. Fe
orbital_total_sum = {}    # e.g. d

# === Read all PDOS files ===
pdos_files = sorted(glob.glob("ppv.pdos_atm#*_wfc#*"))

for file in pdos_files:
    try:
        data = np.loadtxt(file, comments="#")
    except Exception as e:
        print(f"[Skipped] Could not read {file}: {e}")
        continue

    if data.shape[1] < 3:
        print(f"[Skipped] Unexpected format in {file}")
        continue

    energy = data[:, 0]
    up = data[:, 1]
    down = data[:, 2]

    # Extract element and orbital from filename
    filename = os.path.basename(file)
    parts = filename.split("_")
    atom_part = parts[1]  # e.g. atm#1(Fe)
    wfc_part = parts[2]   # e.g. wfc#5(d)

    if "(" not in atom_part or ")" not in atom_part:
        continue
    if "(" not in wfc_part or ")" not in wfc_part:
        continue

    element = atom_part.split("(")[-1].strip(")")
    orbital = wfc_part.split("(")[-1].strip(")")

    if orbital not in ["s", "p", "d"]:
        continue

    # 1. Accumulate element + orbital (e.g. Fe_d)
    key_eo = f"{element}_{orbital}"
    if key_eo not in element_orbital_sum:
        element_orbital_sum[key_eo] = {
            "energy": energy,
            "up": up.copy(),
            "down": down.copy()
        }
    else:
        element_orbital_sum[key_eo]["up"] += up
        element_orbital_sum[key_eo]["down"] += down

    # 2. Accumulate element total (e.g. Fe)
    if element not in element_total_sum:
        element_total_sum[element] = {
            "energy": energy,
            "up": up.copy(),
            "down": down.copy()
        }
    else:
        element_total_sum[element]["up"] += up
        element_total_sum[element]["down"] += down

    # 3. Accumulate orbital total (e.g. d)
    if orbital not in orbital_total_sum:
        orbital_total_sum[orbital] = {
            "energy": energy,
            "up": up.copy(),
            "down": down.copy()
        }
    else:
        orbital_total_sum[orbital]["up"] += up
        orbital_total_sum[orbital]["down"] += down

# === Write element+orbital PDOS (Fe_d_total_pdos.dat) ===
for key, val in element_orbital_sum.items():
    energy, up, down = val["energy"], val["up"], val["down"]
    output = np.column_stack((energy, up, down))
    filename = os.path.join(output_dir, f"{key}_total_pdos.dat")
    header = f"# Fermi energy (eV): {fermi_energy:.6f}\n# Energy (eV)  PDOS_up  PDOS_down"
    np.savetxt(filename, output, fmt="%.6f", header=header, comments="")
    print(f"Saved: {filename}")

# === Write element total PDOS (Fe_total_pdos.dat) ===
for element, val in element_total_sum.items():
    energy, up, down = val["energy"], val["up"], val["down"]
    output = np.column_stack((energy, up, down))
    filename = os.path.join(output_dir, f"{element}_total_pdos.dat")
    header = f"# Fermi energy (eV): {fermi_energy:.6f}\n# Energy (eV)  PDOS_up  PDOS_down"
    np.savetxt(filename, output, fmt="%.6f", header=header, comments="")
    print(f"Saved: {filename}")

# === Write orbital total PDOS (p_total_pdos.dat, etc.) ===
for orbital, val in orbital_total_sum.items():
    energy, up, down = val["energy"], val["up"], val["down"]
    output = np.column_stack((energy, up, down))
    filename = os.path.join(output_dir, f"{orbital}_total_pdos.dat")
    header = f"# Fermi energy (eV): {fermi_energy:.6f}\n# Energy (eV)  PDOS_up  PDOS_down"
    np.savetxt(filename, output, fmt="%.6f", header=header, comments="")
    print(f"Saved: {filename}")

# === Process total DOS (ppv.pdos_tot) ===
dos_file = "ppv.pdos_tot"
if os.path.isfile(dos_file):
    try:
        total_dos = np.loadtxt(dos_file, comments="#")
        energy = total_dos[:, 0]
        up = total_dos[:, 1]
        down = total_dos[:, 2]
        output = np.column_stack((energy, up, down))
        filename = os.path.join(output_dir, "total_DOS.dat")
        header = f"# Fermi energy (eV): {fermi_energy:.6f}\n# Energy (eV)  DOS_up  DOS_down"
        np.savetxt(filename, output, fmt="%.3f", header=header, comments="")
        print(f"Saved: {filename}")
    except Exception as e:
        print(f"[Error] Could not read total DOS file: {e}")
else:
    print("[Warning] total DOS file not found: ppv.pdos_tot")

print("✅ All processing completed.")
