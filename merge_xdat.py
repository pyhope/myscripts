#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path

def parse_xdatcar_header(filename):
    """
    Parse the XDATCAR header and return:
        comment_line
        scale_factor_line
        lattice_lines (3 lines as strings)
        element_line
        count_line
        n_atoms
        header_end_index  (line index where the first configuration starts)
        lattice_values    (3x3 float matrix)
        element_symbols   (list of str)
        element_counts    (list of int)
    """
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) < 8:
        raise ValueError(f"{filename}: file is too short to be a valid XDATCAR.")

    comment_line = lines[0]
    scale_factor_line = lines[1]
    lattice_lines = lines[2:5]
    element_line = lines[5]
    count_line = lines[6]

    try:
        lattice_values = []
        for i in range(3):
            parts = lattice_lines[i].split()
            if len(parts) < 3:
                raise ValueError
            lattice_values.append([float(parts[0]), float(parts[1]), float(parts[2])])

        element_symbols = element_line.split()
        element_counts = [int(x) for x in count_line.split()]
    except Exception as exc:
        raise ValueError(f"{filename}: failed to parse header.") from exc

    if len(element_symbols) == 0:
        raise ValueError(f"{filename}: no element symbols found.")
    if len(element_symbols) != len(element_counts):
        raise ValueError(
            f"{filename}: number of element symbols does not match number of counts."
        )

    n_atoms = sum(element_counts)
    header_end_index = 7

    return (
        comment_line,
        scale_factor_line,
        lattice_lines,
        element_line,
        count_line,
        n_atoms,
        header_end_index,
        lattice_values,
        element_symbols,
        element_counts,
        lines,
    )


def floats_close(a, b, tol=1e-6):
    return abs(a - b) <= tol


def compare_headers(ref, cur, ref_name, cur_name, tol=1e-6):
    """
    Compare scale factor, lattice vectors, element symbols/order, and counts.
    """
    ref_scale = float(ref["scale_factor_line"].split()[0])
    cur_scale = float(cur["scale_factor_line"].split()[0])

    if not floats_close(ref_scale, cur_scale, tol):
        raise ValueError(
            f"Scale factor mismatch:\n"
            f"  {ref_name}: {ref_scale}\n"
            f"  {cur_name}: {cur_scale}"
        )

    for i in range(3):
        for j in range(3):
            a = ref["lattice_values"][i][j]
            b = cur["lattice_values"][i][j]
            if not floats_close(a, b, tol):
                raise ValueError(
                    f"Lattice mismatch at vector {i+1}, component {j+1}:\n"
                    f"  {ref_name}: {a}\n"
                    f"  {cur_name}: {b}"
                )

    if ref["element_symbols"] != cur["element_symbols"]:
        raise ValueError(
            f"Element symbols/order mismatch:\n"
            f"  {ref_name}: {' '.join(ref['element_symbols'])}\n"
            f"  {cur_name}: {' '.join(cur['element_symbols'])}"
        )

    if ref["element_counts"] != cur["element_counts"]:
        raise ValueError(
            f"Element counts mismatch:\n"
            f"  {ref_name}: {' '.join(map(str, ref['element_counts']))}\n"
            f"  {cur_name}: {' '.join(map(str, cur['element_counts']))}"
        )


def extract_configurations(lines, n_atoms, start_index, filename):
    """
    Extract all configurations from XDATCAR.
    Return a list of configurations, each item is:
        (original_header_line, atom_lines)
    where atom_lines has length n_atoms.
    """
    configs = []
    i = start_index
    n_lines = len(lines)

    while i < n_lines:
        line = lines[i].strip()

        if line == "":
            i += 1
            continue

        if not line.startswith("Direct configuration="):
            raise ValueError(
                f"{filename}: unexpected line at {i+1}: {lines[i].rstrip()}\n"
                f"Expected a line starting with 'Direct configuration='."
            )

        if i + n_atoms >= n_lines:
            raise ValueError(
                f"{filename}: incomplete configuration starting at line {i+1}."
            )

        atom_lines = lines[i + 1:i + 1 + n_atoms]
        if len(atom_lines) != n_atoms:
            raise ValueError(
                f"{filename}: configuration at line {i+1} has insufficient atom lines."
            )

        configs.append((lines[i], atom_lines))
        i += 1 + n_atoms

    return configs


def build_header_dict(parsed):
    return {
        "comment_line": parsed[0],
        "scale_factor_line": parsed[1],
        "lattice_lines": parsed[2],
        "element_line": parsed[3],
        "count_line": parsed[4],
        "n_atoms": parsed[5],
        "header_end_index": parsed[6],
        "lattice_values": parsed[7],
        "element_symbols": parsed[8],
        "element_counts": parsed[9],
        "lines": parsed[10],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple VASP XDATCAR files after checking lattice, element order, and counts."
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input XDATCAR files to merge, in the desired order."
    )
    parser.add_argument(
        "-o", "--output",
        default="XDATCAR_merged",
        help="Output merged XDATCAR filename (default: XDATCAR_merged)."
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-6,
        help="Tolerance for comparing lattice parameters and scale factor (default: 1e-6)."
    )

    args = parser.parse_args()
    tol = args.tol

    if len(args.input_files) < 2:
        print("Warning: only one input file provided. The script will still rewrite it with renumbered configurations.", file=sys.stderr)

    input_paths = [Path(f) for f in args.input_files]
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(f"Input file not found: {path}")

    parsed_list = []
    for path in input_paths:
        parsed = parse_xdatcar_header(path)
        parsed_list.append((str(path), build_header_dict(parsed)))

    ref_name, ref = parsed_list[0]

    for cur_name, cur in parsed_list[1:]:
        compare_headers(ref, cur, ref_name, cur_name, tol=tol)

    all_configs = []
    for fname, data in parsed_list:
        configs = extract_configurations(
            lines=data["lines"],
            n_atoms=data["n_atoms"],
            start_index=data["header_end_index"],
            filename=fname
        )
        if len(configs) == 0:
            raise ValueError(f"{fname}: no configurations found.")
        all_configs.extend(configs)

    with open(args.output, "w", encoding="utf-8") as fout:
        fout.write(ref["comment_line"])
        fout.write(ref["scale_factor_line"])
        for line in ref["lattice_lines"]:
            fout.write(line)
        fout.write(ref["element_line"])
        fout.write(ref["count_line"])

        for idx, (_, atom_lines) in enumerate(all_configs, start=1):
            fout.write(f"Direct configuration={idx:6d}\n")
            for line in atom_lines:
                fout.write(line)

    print(f"Merged {len(input_paths)} XDATCAR files into: {args.output}")
    print(f"Total configurations written: {len(all_configs)}")


if __name__ == "__main__":
    main()
