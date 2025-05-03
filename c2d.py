#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("--input_file", "-i", type=str, default="POSCAR",  help="")
parser.add_argument("--output_file", "-o", type=str, default="conf.vasp",  help="")
args = parser.parse_args()

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    system_name = lines[0].strip()
    scale = float(lines[1].strip())

    lattice_vectors = np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
    lattice_vectors *= scale  # apply scaling factor

    element_names = lines[5].split()
    element_counts = list(map(int, lines[6].split()))

    coord_type = lines[7].strip().lower()
    if coord_type[0] not in ['c', 'd']:
        raise ValueError("Coordinate system must be either Cartesian or Direct.")

    num_atoms = sum(element_counts)
    coords = np.array([list(map(float, lines[i].split())) for i in range(8, 8 + num_atoms)])

    return {
        'system': system_name,
        'lattice': lattice_vectors,
        'elements': element_names,
        'counts': element_counts,
        'coord_type': coord_type,
        'coordinates': coords
    }

def cartesian_to_fractional(cart_coords, lattice):
    inv_lattice = np.linalg.inv(lattice.T)
    return np.dot(cart_coords, inv_lattice)

def write_poscar(data, frac_coords, output_filename='POSCAR_direct'):
    with open(output_filename, 'w') as f:
        f.write(f"{data['system']}\n")
        f.write("1.000000\n")
        for vec in data['lattice']:
            f.write(f"  {vec[0]:.9f}  {vec[1]:.9f}  {vec[2]:.9f}\n")
        f.write("  " + "  ".join(data['elements']) + "\n")
        f.write("  " + "  ".join(map(str, data['counts'])) + "\n")
        f.write("Direct\n")
        for coord in frac_coords:
            f.write(f"  {coord[0]:.9f}  {coord[1]:.9f}  {coord[2]:.9f}\n")

poscar_file = args.input_file
data = read_poscar(poscar_file)
if data['coord_type'] == 'cartesian':
    frac_coords = cartesian_to_fractional(data['coordinates'], data['lattice'])
    write_poscar(data, frac_coords, output_filename=args.output_file)
    print("Conversion completed. Output written to 'POSCAR_direct'.")
else:
    print("POSCAR is already in Direct format. No conversion needed.")
