#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial.distance import cdist
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--input_file", "-i", type=str, default="conf.qe")
parser.add_argument("--output_file", "-o", type=str, default="sorted.qe")
parser.add_argument("--v_output_file", "-v", type=str, default="V.txt")
parser.add_argument("--cutoff", "-c", type=float, default=4.0)
parser.add_argument("--k", type=float, default=-0.6373)
parser.add_argument("--m", type=float, default=3.0547)

args = parser.parse_args()

def V_func(d, k=args.k, m=args.m):
    return k * d + m

def read_qe_structure(file_path):
    """ Read Quantum Espresso structure file and extract cell parameters and atomic positions """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    cell_start = next(i for i, line in enumerate(lines) if line.startswith("CELL_PARAMETERS"))
    atomic_start = next(i for i, line in enumerate(lines) if line.startswith("ATOMIC_POSITIONS"))
    
    cell = np.array([list(map(float, lines[cell_start + i + 1].split())) for i in range(3)])
    atoms = [line.split() for line in lines[atomic_start + 1:]]
    
    return cell, atoms

def fractional_to_cartesian(positions, cell):
    """ Convert atomic positions from fractional to Cartesian coordinates """
    return np.dot(np.array(positions), cell)

def compute_min_distances(atoms, cell):
    """Compute the minimum distance of each O to any Fe considering PBC."""
    fractional_positions = np.array([list(map(float, atom[1:])) for atom in atoms])
    cartesian_positions = fractional_to_cartesian(fractional_positions, cell)

    fe_positions = cartesian_positions[[i for i, atom in enumerate(atoms) if atom[0] == 'Fe']]
    o_indices = [i for i, atom in enumerate(atoms) if atom[0] == 'O']

    if len(fe_positions) == 0 or len(o_indices) == 0:
        return {i: float('inf') for i in o_indices}

    min_distances = {}
    print("Calculating distances for O atoms...")

    for i, o_idx in enumerate(o_indices):
        o_pos = cartesian_positions[o_idx]
        min_dist = float('inf')

        for shift in np.array([[i, j, k] for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1]]):
            shifted_fe_positions = fe_positions + np.dot(shift, cell)
            distances = cdist([o_pos], shifted_fe_positions)
            min_dist = min(min_dist, np.min(distances))

        min_distances[o_idx] = min_dist
        print(f"O atom {o_idx+1}: min distance to Fe = {min_dist:.4f}")

    sorted_o_atoms = sorted(min_distances.items(), key=lambda x: x[1])
    print("Top 10 closest O atoms:")
    for idx, (o_idx, dist) in enumerate(sorted_o_atoms):
        print(f"{idx+1}: O atom index {o_idx+1}, Distance {dist:.4f}")

    return min_distances

def sort_atoms(atoms, min_distances):
    """Sort atoms: Fe first, then O sorted by distance, then other elements in original order."""
    fe_atoms = [atom for atom in atoms if atom[0] == 'Fe']
    o_atoms = sorted((atom for i, atom in enumerate(atoms) if atom[0] == 'O'),
                     key=lambda x: min_distances[atoms.index(x)])
    other_atoms = [atom for atom in atoms if atom[0] not in {'Fe', 'O'}]

    return fe_atoms + o_atoms + other_atoms

def write_qe_structure(output_path, cell, atoms):
    """ Write the sorted structure back to a Quantum Espresso format file """
    with open(output_path, 'w') as f:
        f.write("CELL_PARAMETERS (angstrom)\n")
        for row in cell:
            f.write("  " + "  ".join(f"{x:.8f}" for x in row) + "\n")
        f.write("\nATOMIC_POSITIONS (crystal)\n")
        for atom in atoms:
            f.write("  " + "  ".join(atom) + "\n")
    print(f"Sorted structure written to {output_path}")

def find_fe_o_distances(atoms, cell, cutoff=args.cutoff):
    """Find Fe-O pairs within cutoff distance and compute V values."""
    fractional_positions = np.array([list(map(float, atom[1:])) for atom in atoms])
    cartesian_positions = fractional_to_cartesian(fractional_positions, cell)
    
    fe_indices = [i for i, atom in enumerate(atoms) if atom[0] == 'Fe']
    o_indices = [i for i, atom in enumerate(atoms) if atom[0] == 'O']
    
    output_lines = []
    print_lines = []
    
    print(f"Found {len(fe_indices)} Fe atoms and {len(o_indices)} O atoms.")
    
    for fe_idx in fe_indices:
        fe_pos = cartesian_positions[fe_idx]
        
        for o_idx in o_indices:
            o_pos = cartesian_positions[o_idx]
            
            min_dist = float('inf')
            for shift in np.array([[i, j, k] for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1]]):
                shifted_o_pos = o_pos + np.dot(shift, cell)
                dist = np.linalg.norm(fe_pos - shifted_o_pos)
                min_dist = min(min_dist, dist)
            
            if min_dist < cutoff:
                fe_number = fe_idx + 1  # 1-based index
                o_number = o_idx + 1    # 1-based index
                v_value = round(V_func(min_dist), 2)
                output_lines.append(f"V Fe-3d O-2p {fe_number} {o_number} {v_value}")
                print_lines.append(f"Fe{fe_number}-O{o_number}: d = {min_dist:.4f}, V = {v_value}")
        
        print(f"Processed Fe atom {fe_number}, found {len(output_lines)} valid Fe-O pairs so far.")
    
    return print_lines, output_lines

input_file = args.input_file
output_file = args.output_file

print(f"Reading structure from {input_file}...")
cell, atoms = read_qe_structure(input_file)
min_distances = compute_min_distances(atoms, cell)
print("Sorting atoms...")
sorted_atoms = sort_atoms(atoms, min_distances)
write_qe_structure(output_file, cell, sorted_atoms)
print("Process completed successfully.")
print("*" * 80)
print("\nFinding Fe-O distances...")
print_lines, output_lines = find_fe_o_distances(sorted_atoms, cell)

print("\nGenerated output:")
for line in print_lines:
    print(line)

with open(args.v_output_file, 'w') as f:
    f.write("\n".join(output_lines))
