#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.dont_write_bytecode = True
import numpy as np
import argparse
import csv
import random
from c2d import read_poscar, fractional_to_cartesian, cartesian_to_fractional

def pbc_distance(vec1, vec2, lattice):
    # Compute periodic distance under PBC using minimum image convention
    delta = vec2 - vec1
    frac = np.dot(delta, np.linalg.inv(lattice.T))
    frac -= np.round(frac)
    cart = np.dot(frac, lattice.T)
    return np.linalg.norm(cart)

def main():
    parser = argparse.ArgumentParser(description="Replace top-n A atoms with smallest d1 and largest d2 to B, output structure in Direct coordinates.")
    parser.add_argument("--input_file", "-i", type=str, default="POSCAR", help="Input POSCAR file")
    parser.add_argument("--output_file", "-o", type=str, default="modified.vasp", help="Output modified structure file")
    parser.add_argument("--element_a", "-a", type=str, default="Si", help="Target atom A (default: Si)")
    parser.add_argument("--element_b", "-b", type=str, default="Fe", help="Neighbor atom B (default: Fe)")
    parser.add_argument("--element_c", "-c", type=str, default="Al", help="Replacement element C (default: Al)")
    parser.add_argument("--element_d", "-d", type=str, default="Mg", help="Atoms after replacement during sorting (default: Mg)")
    parser.add_argument("--top_n", "-n", type=int, default=None, help="Number of A atoms to replace (default: half of B atoms)")
    parser.add_argument("--save_csv", "-s", action="store_true", default=False, help="Save replacement pairs to CSV file")
    parser.add_argument("--csv_filename", "-csv", type=str, default="replaced_pairs.csv", help="CSV file for replacement pairs (default: replaced_pairs.csv)")
    parser.add_argument("--disable_high_d2", "-dhd", action="store_true", default=False, help="Disable high d2 filtering")
    parser.add_argument("--random_seed", "-rs", type=int, default=None, help="Random seed for reproducibility (default: None)")
    args = parser.parse_args()

    # Load structure data
    data = read_poscar(args.input_file)
    lattice = data["lattice"]
    coord_type = data["coord_type"].lower()

    if coord_type == "cartesian":
        coords_cart = data["coordinates"]
    elif coord_type == "direct":
        print("Converting Direct coordinates to Cartesian...")
        coords_cart = fractional_to_cartesian(data["coordinates"], lattice)
    else:
        raise ValueError("Unknown coordinate type: must be Cartesian or Direct.")

    elements = data["elements"]
    counts = data["counts"]
    comments = data["comments"]

    # Build element list by atom index
    element_list = []
    for elem, count in zip(elements, counts):
        element_list.extend([elem] * count)

    indices = list(range(len(element_list)))
    indices_a = [i for i, e in enumerate(element_list) if e == args.element_a]
    indices_b = [i for i, e in enumerate(element_list) if e == args.element_b]
    indices_d = [i for i, e in enumerate(element_list) if e == args.element_d]

    if args.top_n is None:
        n_replace = len(indices_b) // 2
    else:
        n_replace = args.top_n

    tol = 1e-3
    all_pairs_info = []
    d1_values = []
    a_d_distances = []

    # Analyze A atoms
    for ia in indices_a:
        a_pos = coords_cart[ia]
        distances = []
        
        for id in indices_d:
            dist = pbc_distance(a_pos, coords_cart[id], lattice)
            a_d_distances.append(dist)

        for ib in indices_b:
            b_pos = coords_cart[ib]
            dist = pbc_distance(a_pos, b_pos, lattice)
            distances.append((dist, ib))

        if not distances:
            continue

        distances.sort()
        d1 = distances[0][0]
        nearest = [d for d in distances if abs(d[0] - d1) < tol]

        if len(nearest) > 1:
            continue  # ambiguous nearest B

        b1_index = nearest[0][1]
        second_neighbors = [d for d in distances if abs(d[0] - d1) >= tol]
        if not second_neighbors:
            continue
        d2 = second_neighbors[0][0]

        all_pairs_info.append((d1, d2, ia, b1_index))
        d1_values.append(d1)

    if not all_pairs_info:
        print("No unique nearest neighbors found.")
        return

    # Get global minimum distance
    d1_min = np.min(a_d_distances)

    # Filter A atoms with d1 ≈ d1_min
    d2_candidates = [(d2, d1, ia, ib) for (d1, d2, ia, ib) in all_pairs_info if abs(d1 - d1_min) < tol]
    if not d2_candidates:
        print(f"No {args.element_a} atoms found with d1 ≈ d1_min ({d1_min:.6f} Å).")
        return

    min_d2 = min(d2 for d2, _, _, _ in d2_candidates)

    if args.disable_high_d2:
        filtered_candidates = d2_candidates
    else:
        # Filter candidates with d2 significantly larger than the minimum
        filtered_candidates = [(d2, d1, ia, ib) for (d2, d1, ia, ib) in d2_candidates if d2 > min_d2 + tol]
        if not filtered_candidates:
            print("No candidates with d2 significantly larger than the minimum.")
            return

    if args.random_seed is not None:
        random.seed(args.random_seed)
    random.shuffle(filtered_candidates)

    selected = []
    used_b1_indices = set()

    for d2, d1, ia, ib in filtered_candidates:
        if ib in used_b1_indices:
            continue  # Skip if B1 already used
        selected.append((d2, d1, ia, ib))
        used_b1_indices.add(ib)
        if len(selected) == n_replace:
            break

    if len(selected) < n_replace:
        print(f"Warning: Only {len(selected)} {args.element_a} atoms found with d1 ≈ d1_min. Replacing these.")
        n_replace = len(selected)

    selected_a_indices = [ia for _, _, ia, _ in selected]
    replaced_b1_indices = [ib for _, _, _, ib in selected]

    # Replace selected A atoms with C
    for ia in selected_a_indices:
        element_list[ia] = args.element_c

    # Write replaced_pairs.csv
    if args.save_csv:
        with open(args.csv_filename, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow([f"{args.element_a} index", "x", "y", "z", f"{args.element_b} index", "x", "y", "z", "d1 (Å)", "d2 (Å)", "Replaced"])

            replaced_rows = []
            non_replaced_rows = []

            for d1, d2, ia, ib in all_pairs_info:
                row = [
                    ia, *coords_cart[ia], ib, *coords_cart[ib], f"{d1:.6f}", f"{d2:.6f}",
                    "Yes" if ia in selected_a_indices else "No"
                ]
                if ia in selected_a_indices:
                    replaced_rows.append(row)
                else:
                    non_replaced_rows.append(row)

            # Write replaced first, then the rest
            for row in replaced_rows + non_replaced_rows:
                writer.writerow(row)
            print(f"Replacement info saved to '{args.csv_filename}'.")

    # Sort indices for output POSCAR
    replaced_b1_set = set(replaced_b1_indices)
    replaced_a_set = set(selected_a_indices)

    sorted_indices = []
    for i in indices_b:
        if i in replaced_b1_set:
            sorted_indices.append(i)
            comments[i] += f" => {args.element_b}3+"
    for i in indices_b:
        if i not in replaced_b1_set:
            sorted_indices.append(i)
            comments[i] += f" => {args.element_b}2+"
    for i in selected_a_indices:
            sorted_indices.append(i)
            comments[i] += f" => {args.element_c}"
    sorted_indices += [i for i in indices_d]
    sorted_indices += [i for i in indices_a if i not in replaced_a_set]
    sorted_indices += [i for i in indices if i not in indices_b and i not in indices_a and i not in indices_d]

    sorted_coords_cart = [coords_cart[i] for i in sorted_indices]
    sorted_coords_direct = cartesian_to_fractional(np.array(sorted_coords_cart), lattice)
    sorted_comments = [comments[i] for i in sorted_indices]
    sorted_elements = [element_list[i] for i in sorted_indices]

    # Collect unique elements and counts in order
    new_elements = []
    new_counts = []

    for elem in [args.element_b, args.element_c, args.element_d, args.element_a]:
        count = sorted_elements.count(elem)
        if count > 0:
            new_elements.append(elem)
            new_counts.append(count)

    for elem in elements:
        if elem not in new_elements:
            count = sorted_elements.count(elem)
            if count > 0:
                new_elements.append(elem)
                new_counts.append(count)

    # Write modified.vasp in Direct coordinates
    with open(args.output_file, 'w') as f:
        f.write(f"{data['system']} - modified\n")
        f.write("1.000000\n")
        for vec in lattice:
            f.write(f"  {vec[0]:.9f}  {vec[1]:.9f}  {vec[2]:.9f}\n")
        f.write("  " + "  ".join(new_elements) + "\n")
        f.write("  " + "  ".join(map(str, new_counts)) + "\n")
        f.write("Direct\n")
        for coord, comment in zip(sorted_coords_direct, sorted_comments):
            f.write(f"  {coord[0]:.9f}  {coord[1]:.9f}  {coord[2]:.9f}  {comment}\n")

    print(f"Replaced {len(selected_a_indices)} '{args.element_a}' atoms with '{args.element_c}'.")
    print(f"New structure written to '{args.output_file}' in Direct coordinates.")

if __name__ == "__main__":
    main()
