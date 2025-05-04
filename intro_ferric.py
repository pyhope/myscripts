#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import csv
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
    parser.add_argument("--top_n", "-n", type=int, default=None, help="Number of A atoms to replace (default: half of B atoms)")
    parser.add_argument("--save_csv", "-s", action="store_true", default=False, help="Save replacement pairs to CSV file")
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

    if args.top_n is None:
        n_replace = len(indices_b) // 2
    else:
        n_replace = args.top_n

    tol = 1e-3
    all_pairs_info = []
    d1_values = []

    # Analyze A atoms
    for ia in indices_a:
        a_pos = coords_cart[ia]
        distances = []

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

    # Get global minimum d1
    d1_min = min(d1_values)

    # Filter A atoms with d1 ≈ d1_min
    d2_candidates = [(d2, d1, ia, ib) for (d1, d2, ia, ib) in all_pairs_info if abs(d1 - d1_min) < tol]
    d2_candidates.sort(reverse=True)
    selected = d2_candidates[:n_replace]

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
        with open("replaced_pairs.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["A_index", "A_x", "A_y", "A_z", "B1_index", "B1_x", "B1_y", "B1_z", "d1", "d2", "Replaced"])
            for d1, d2, ia, ib in all_pairs_info:
                replaced_flag = "Yes" if ia in selected_a_indices else "No"
                writer.writerow([
                    ia, *coords_cart[ia], ib, *coords_cart[ib], f"{d1:.6f}", f"{d2:.6f}", replaced_flag
                ])

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
    sorted_indices += [i for i in indices_a if i not in replaced_a_set]
    sorted_indices += [i for i in indices if i not in indices_b and i not in indices_a]

    sorted_coords_cart = [coords_cart[i] for i in sorted_indices]
    sorted_coords_direct = cartesian_to_fractional(np.array(sorted_coords_cart), lattice)
    sorted_comments = [comments[i] for i in sorted_indices]
    sorted_elements = [element_list[i] for i in sorted_indices]

    # Collect unique elements and counts in order
    new_elements = []
    new_counts = []

    for elem in [args.element_b, args.element_c, args.element_a]:
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
    print("Replacement info saved to 'replaced_pairs.csv'.")

if __name__ == "__main__":
    main()
