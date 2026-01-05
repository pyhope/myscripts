#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys

def move_atom_pbc_min_image(
    poscar_in="relaxed.vasp",
    poscar_out="out.vasp",
    atom_index=103,     # 1-based
    ref_index=2,        # 1-based
    shift=0.3          # Angstrom
):
    with open(poscar_in, "r") as f:
        lines = f.readlines()

    # --- Parse POSCAR (tailored to your example: no Selective Dynamics, Cartesian coords) ---
    # Layout:
    # 0: comment
    # 1: scale
    # 2-4: lattice vectors (Cartesian)
    # 5: element symbols (unused here)
    # 6: element counts -> total atoms
    # 7: "Cartesian" (or "Direct")
    # 8..: coordinates
    scale = float(lines[1].strip())
    A = np.array([
        list(map(float, lines[2].split())),
        list(map(float, lines[3].split())),
        list(map(float, lines[4].split()))
    ]) * scale  # apply global scale

    # total number of atoms
    counts = list(map(int, lines[6].split()))
    n_atoms = sum(counts)

    coord_mode = lines[7].strip().lower()
    if not coord_mode.startswith("cart"):
        raise ValueError("This script expects Cartesian coordinates on line 8.")

    # read coordinates (first n_atoms lines after line 7)
    coord_start = 8
    coord_end = coord_start + n_atoms
    coords = []
    for i in range(coord_start, coord_end):
        parts = lines[i].split()
        coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
    coords = np.array(coords)  # Cartesian (Å)

    # --- Minimum-image displacement under PBC ---
    # Convert to fractional, wrap delta to [-0.5, 0.5), then back to Cartesian
    A_inv = np.linalg.inv(A)

    r_i = coords[atom_index - 1]           # target atom (to move)
    r_j = coords[ref_index - 1]            # reference atom
    f_i = A_inv @ r_i
    f_j = A_inv @ r_j

    df = f_j - f_i
    df_wrapped = df - np.round(df)         # wrap each component into [-0.5, 0.5)

    dcart = A @ df_wrapped                 # shortest Cartesian displacement
    norm = np.linalg.norm(dcart)
    if norm < 1e-12:
        raise ValueError("Minimum-image displacement is zero; cannot define a direction.")

    direction = dcart / norm
    coords[atom_index - 1] = r_i + shift * direction

    # --- Write out POSCAR_new (preserve header; overwrite coord block) ---
    out_lines = lines[:coord_start]
    for c in coords:
        out_lines.append(f"  {c[0]:.9f}  {c[1]:.9f}  {c[2]:.9f}\n")
    # If the input had exactly n_atoms coordinate lines and nothing after, we’re done.
    # If there is extra content after coordinates, you could append lines[coord_end:].

    with open(poscar_out, "w") as f:
        f.writelines(out_lines)

    print(f"Atom {atom_index} moved by {shift} Å towards the nearest periodic image of atom {ref_index}.")
    print(f"Output written to {poscar_out}")

if __name__ == "__main__":
    move_atom_pbc_min_image(
        poscar_out="ini.vasp",
        atom_index=int(sys.argv[1]),
        ref_index=1,
    )
    move_atom_pbc_min_image(
        poscar_out="fin.vasp",
        atom_index=int(sys.argv[2]),
        ref_index=2,
    )
