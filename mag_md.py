#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--Fe_idx", "-i", type=int, default=160,  help="Index of the Fe atom")
parser.add_argument("--directory_name", "-d", type=str, default="./normal",  help="Directory name")

args = parser.parse_args()

Fe_idx = args.Fe_idx
dir_name = args.directory_name

def process_outcar():
    output_lines = []
    
    # Read all lines from the OUTCAR file
    with open(str(dir_name) + "/OUTCAR", "r") as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        # Strip whitespace from both ends and split into words
        tokens = lines[i].strip().split()
        # Check if there are enough words and if the second and third words are "Ionic" and "step"
        if len(tokens) >= 4 and tokens[1] == "Ionic" and tokens[2] == "step":
            # The fourth word represents the step number
            step = tokens[3]
            
            # Look for the line containing magnetization (x)
            magnetic_moment = None
            j = i + 1
            found_magnet_line = False
            while j < len(lines):
                tokens_j = lines[j].strip().split()
                # Check if the first two words are "magnetization" and "(x)"
                if len(tokens_j) >= 2 and tokens_j[0] == "magnetization" and tokens_j[1] == "(x)":
                    found_magnet_line = True
                    break
                j += 1
            
            # If the magnetization (x) line was found, search from that line for the first line starting with "397"
            if found_magnet_line:
                k = j
                while k < len(lines):
                    tokens_k = lines[k].strip().split()
                    if len(tokens_k) > 2:
                        if tokens_k[0] == str(Fe_idx):
                            magnetic_moment = tokens_k[-1]
                            break
                    k += 1
            
            # Look for the line starting with "FREE ENERGIE" and record the 5th word as the energy
            energy = None
            m = i + 1
            while m < len(lines):
                tokens_m = lines[m].strip().split()
                if len(tokens_m) >= 5 and tokens_m[0] == "FREE" and tokens_m[1] == "ENERGIE":
                    toten = lines[m+2].strip().split()[4]
                    energy = lines[m+4].strip().split()[3]
                    break
                m += 1
            
            # If all values are found, append the step number, toten, energy, and magnetic moment to the output list
            if step is not None and energy is not None and magnetic_moment is not None:
                output_lines.append([step, toten, energy, magnetic_moment])
        
        i += 1

    # Write the results to the output.txt file in the current directory
    for line in output_lines:
        # if 150 <= int(line[0]) <=195 and int(line[0]) % 5 == 0:
            print(" ".join(line))

if __name__ == "__main__":
    process_outcar()
