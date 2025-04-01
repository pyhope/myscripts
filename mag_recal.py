#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--Fe_idx", "-i", type=int, default=160,  help="Index of the Fe atom")
parser.add_argument("--subdirectory_name", "-d", type=str, default="./8",  help="Subdirectory name")
parser.add_argument("--directory_num", "-n", type=int, default=10,  help="Directory number")

args = parser.parse_args()

Fe_idx = args.Fe_idx
dir_name = args.subdirectory_name
dir_num = args.directory_num

def process_outcar_in_dir(directory):
    """
    Process the OUTCAR file in the specified directory:
      - Use the directory name as the ionic step value.
      - Find the line containing "magnetization (x)", then, starting from that line, search downward until the first line starting with "395" is encountered,
        and record the last word of that line as the magnetic moment.
      - Traverse the entire file to find the last occurrence of the line containing "TOTEN", and use the second-to-last word of that line as the energy.
    If both energy and magnetic moment are found, return a formatted list [step, energy, energy2, magnetic_moment]; otherwise, return None.
    """
    outcar_path = os.path.join(directory, dir_name + "/OUTCAR")
    try:
        with open(outcar_path, "r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Unable to open {outcar_path}: {e}")
        return None

    magnetic_moment = None
    energy = None

    # Find the line containing "magnetization (x)"
    magnet_line_index = None
    for idx, line in enumerate(lines):
        tokens = line.strip().split()
        if len(tokens) >= 2 and tokens[0] == "magnetization" and tokens[1] == "(x)":
            magnet_line_index = idx
            break

    # Starting from the "magnetization (x)" line, search for the first line that begins with "397" and record its last word as the magnetic moment
    if magnet_line_index is not None:
        for line in lines[magnet_line_index:]:
            tokens = line.strip().split()
            if tokens and tokens[0] == str(Fe_idx):
                magnetic_moment = tokens[-1]
                break

    # Traverse the entire OUTCAR file and record the last occurrence of a line containing "TOTEN" along with its tokens
    last_toten_index = None
    last_toten_tokens = None
    for idx, line in enumerate(lines):
        if "TOTEN" in line:
            tokens = line.strip().split()
            if "TOTEN" in tokens:
                last_toten_index = idx
                last_toten_tokens = tokens

    if last_toten_tokens is not None and len(last_toten_tokens) >= 2:
        energy = last_toten_tokens[-2]
    else:
        print(f"No valid TOTEN line found in {outcar_path} to extract energy.")

    # From two lines after the last TOTEN line, extract energy2 (the 5th token)
    if last_toten_index is not None and (last_toten_index + 2) < len(lines):
        tokens_energy2 = lines[last_toten_index + 2].strip().split()
        if len(tokens_energy2) >= 5:
            energy2 = tokens_energy2[3]
        else:
            print(f"Not enough tokens in the second line after the TOTEN line in {outcar_path} to extract energy2.")
    else:
        print(f"There are not enough lines after the TOTEN line in {outcar_path} to extract energy2.")

    # Use the directory name as the ionic step value
    step = os.path.basename(directory)

    if energy is not None and energy2 is not None and magnetic_moment is not None:
        return [step, energy, energy2, magnetic_moment]
    else:
        print(f"Incomplete data in directory {directory}")
        return None

def process_all_outcars():
    """
    Traverse subdirectories named 1 through 12 in the current directory.
    Each subdirectory should contain an OUTCAR file.
    After processing, print all records in the format "ionic_step energy energy2 magnetic_moment".
    """
    output_lines = []
    for step in range(1, dir_num + 1):
        dir_name = str(step)
        if os.path.isdir(dir_name):
            result = process_outcar_in_dir(dir_name)
            if result:
                output_lines.append(result)
        else:
            print(f"Directory {dir_name} does not exist.")
    
    for line in output_lines:
        print(" ".join(line))

if __name__ == "__main__":
    process_all_outcars()
