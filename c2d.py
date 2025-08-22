#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = [ln.rstrip("\n") for ln in f]

    if len(lines) < 8:
        raise ValueError("POSCAR too short.")

    system_name = lines[0].strip()
    scale = float(lines[1].strip())

    # lattice as rows: a, b, c
    lattice = np.array([list(map(float, lines[2].split())),
                        list(map(float, lines[3].split())),
                        list(map(float, lines[4].split()))], dtype=float)
    lattice *= scale

    # elements and counts
    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))
    num_atoms = sum(counts)

    # optional "Selective dynamics" line
    line7 = lines[7].strip()
    sel_dyn = False
    coord_type_line_idx = 7
    if line7.lower().startswith("selective"):
        sel_dyn = True
        coord_type_line_idx = 8

    coord_type = lines[coord_type_line_idx].strip().lower()
    if coord_type not in ("direct", "cartesian"):
        raise ValueError(f"Unknown coordinate type: {lines[coord_type_line_idx]!r}")

    start = coord_type_line_idx + 1
    coord_lines = lines[start:start + num_atoms]
    if len(coord_lines) < num_atoms:
        raise ValueError("Not enough coordinate lines.")

    coords = []
    flags = []
    comments = []
    for ln in coord_lines:
        # split off comment starting with '!'
        if '!' in ln:
            coord_part, comment_part = ln.split('!', 1)
            comment = '! ' + comment_part.strip()
        else:
            coord_part = ln
            comment = ''
        parts = coord_part.split()
        if len(parts) < 3:
            raise ValueError(f"Malformed coordinate line: {ln!r}")
        x, y, z = map(float, parts[:3])
        coords.append([x, y, z])

        # selective dynamics flags if present
        if sel_dyn and len(parts) >= 6:
            flags.append(parts[3:6])
        elif sel_dyn:
            # fill default T T T if missing
            flags.append(["T", "T", "T"])
        else:
            flags.append(None)
        comments.append(comment)

    return {
        "system": system_name,
        "lattice": lattice,          # rows = a,b,c
        "elements": elements,
        "counts": counts,
        "coord_type": coord_type,    # 'direct' or 'cartesian'
        "selective": sel_dyn,
        "coordinates": np.array(coords, dtype=float),
        "flags": flags,
        "comments": comments,
    }

def cartesian_to_fractional(cart_coords, lattice):
    """
    cart = frac @ lattice  =>  frac = cart @ inv(lattice)
    """
    inv_lat = np.linalg.inv(lattice)
    return cart_coords @ inv_lat

def fractional_to_cartesian(frac_coords, lattice):
    """
    cart = frac @ lattice
    """
    return frac_coords @ lattice

def write_poscar(data, new_coords, new_type, output_filename, wrap=False):
    with open(output_filename, 'w') as f:
        f.write(f"{data['system']}\n")
        f.write("1.000000\n")
        L = data['lattice']
        for i in range(3):
            f.write(f"  {L[i,0]:.9f}  {L[i,1]:.9f}  {L[i,2]:.9f}\n")
        f.write("  " + "  ".join(data['elements']) + "\n")
        f.write("  " + " ".join(map(str, data['counts'])) + "\n")

        if data['selective']:
            f.write("Selective dynamics\n")
        f.write(f"{new_type.capitalize()}\n")

        for i, (coord, comment) in enumerate(zip(new_coords, data['comments'])):
            c = coord.copy()
            if wrap and new_type.lower() == "direct":
                c = c - np.floor(c)   # wrap to [0,1)
            if data['selective'] and data['flags'][i] is not None:
                flag_str = "  " + "  ".join(data['flags'][i])
            else:
                flag_str = ""
            f.write(f"  {c[0]:.9f}  {c[1]:.9f}  {c[2]:.9f}{flag_str}")
            if comment:
                f.write(f"  {comment}")
            f.write("\n")

def main():
    parser = argparse.ArgumentParser(
        description='Convert between Cartesian and Direct coordinates in a POSCAR-like file.'
    )
    parser.add_argument("-i", "--input_file", type=str, default="POSCAR", help="Input POSCAR file")
    parser.add_argument("-o", "--output_file", type=str, default="conf.vasp", help="Output file")
    parser.add_argument("--wrap", action="store_false", default=True, help="When output is Direct, do not wrap fractional coords into [0,1)")
    args = parser.parse_args()

    try:
        data = read_poscar(args.input_file)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    if data['coord_type'] == 'cartesian':
        converted = cartesian_to_fractional(data['coordinates'], data['lattice'])
        new_type = 'Direct'
    elif data['coord_type'] == 'direct':
        converted = fractional_to_cartesian(data['coordinates'], data['lattice'])
        new_type = 'Cartesian'
    else:
        raise ValueError("Unsupported coordinate type: must be 'Cartesian' or 'Direct'.")

    write_poscar(data, converted, new_type, args.output_file, wrap=args.wrap)
    print(f"[OK] Wrote '{args.output_file}' in {new_type} format.")

if __name__ == '__main__':
    main()
