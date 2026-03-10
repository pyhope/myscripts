#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import argparse


def replace_fe_with_element(
    input_file,
    output_file,
    n_replace=8,
    seed=None,
    keep_k=False,
    target_element="Si",
):
    if seed is not None:
        random.seed(seed)

    with open(input_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError("Input file is too short to be a valid extended XYZ file.")

    # Header
    natoms_line = lines[0].rstrip("\n")
    comment_line = lines[1].rstrip("\n")

    # Remove ':k:R:1' from header only if k is not kept
    if not keep_k:
        comment_line = comment_line.replace(":k:R:1", "")

    atom_lines = lines[2:]

    parsed_atoms = []
    candidate_indices = []

    for i, line in enumerate(atom_lines):
        stripped = line.strip()
        if not stripped:
            continue

        parts = stripped.split()

        # Expect at least: species x y z k
        if len(parts) < 5:
            raise ValueError(f"Line {i+3} does not have enough columns:\n{line}")

        species = parts[0]
        xyz = parts[1:4]
        k_value = parts[-1]

        parsed_atoms.append({
            "species": species,
            "xyz": xyz,
            "k": k_value,
        })

        if species == "Fe" and k_value == "2":
            candidate_indices.append(len(parsed_atoms) - 1)

    if len(candidate_indices) < n_replace:
        raise ValueError(
            f"Only {len(candidate_indices)} Fe atoms with k=2 were found, "
            f"but {n_replace} replacements were requested."
        )

    # Randomly choose Fe(k=2) atoms to replace
    selected = random.sample(candidate_indices, n_replace)

    for idx in selected:
        parsed_atoms[idx]["species"] = target_element

    # Move all replaced atoms (target element) to the end
    non_target_atoms = [atom for atom in parsed_atoms if atom["species"] != target_element]
    target_atoms = [atom for atom in parsed_atoms if atom["species"] == target_element]
    final_atoms = non_target_atoms + target_atoms

    # Write output
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"{len(final_atoms)}\n")
        f.write(comment_line + "\n")
        for atom in final_atoms:
            if keep_k:
                f.write(f"{atom['species']} {' '.join(atom['xyz'])} {atom['k']}\n")
            else:
                f.write(f"{atom['species']} {' '.join(atom['xyz'])}\n")


parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", type=str, default="prox.xyz", help="input .xyz file")
parser.add_argument("--output_file", "-o", type=str, default="output.xyz", help="output .xyz file")
parser.add_argument("--n_replace", "-n", type=int, default=8, help="number of Fe atoms to replace")
parser.add_argument("--seed", type=int, default=None, help="random seed for reproducibility (optional)")
parser.add_argument(
    "--keep_k", "-k",
    action="store_true",
    help="keep the last k column in output and preserve ':k:R:1' in the header",
)
parser.add_argument(
    "--element", "-e",
    type=str,
    default="Si",
    help="target element used to replace Fe (default: Si)",
)

args = parser.parse_args()

replace_fe_with_element(
    input_file=args.input_file,
    output_file=args.output_file,
    n_replace=args.n_replace,
    seed=args.seed,
    keep_k=args.keep_k,
    target_element=args.element,
)

print(f"Done. Output written to: {args.output_file}")
