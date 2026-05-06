#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sort atoms in an xyz file by species order."
    )
    parser.add_argument(
        "-i", "--input",
        default="tmp.xyz",
        help="Input xyz file name (default: tmp.xyz)"
    )
    parser.add_argument(
        "-o", "--output",
        default="tmp.sorted.xyz",
        help="Output xyz file name (default: tmp.sorted.xyz)"
    )
    parser.add_argument(
        "--order",
        nargs="+",
        default=["Mg", "O", "Fe", "W"],
        help="Desired species order, e.g. --order Mg O Fe W"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)
    species_order = args.order

    if len(species_order) != len(set(species_order)):
        raise ValueError(f"Duplicate species found in --order: {species_order}")

    species_rank = {species: i for i, species in enumerate(species_order)}

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    with input_file.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise ValueError("Invalid xyz file: file is too short.")

    natoms = int(lines[0].strip())
    header = lines[1]
    atom_lines = lines[2:]

    if len(atom_lines) != natoms:
        raise ValueError(
            f"Atom count mismatch: first line says {natoms}, but found {len(atom_lines)} atom lines."
        )

    parsed_atoms = []
    for lineno, line in enumerate(atom_lines, start=3):
        parts = line.split()
        if len(parts) < 5:
            raise ValueError(f"Invalid atom line at line {lineno}: {line.strip()}")

        species = parts[0]
        if species not in species_rank:
            raise ValueError(
                f"Unknown species '{species}' at line {lineno}. "
                f"It is not included in --order {species_order}"
            )

        parsed_atoms.append((species_rank[species], line))

    # Stable sort: preserve original order within each species
    parsed_atoms.sort(key=lambda item: item[0])
    sorted_atom_lines = [line for _, line in parsed_atoms]

    with output_file.open("w", encoding="utf-8") as f:
        f.write(f"{natoms}\n")
        f.write(header)
        f.writelines(sorted_atom_lines)

    print(f"Species order: {' '.join(species_order)}")
    print(f"Sorted xyz file written to: {output_file}")


if __name__ == "__main__":
    main()
