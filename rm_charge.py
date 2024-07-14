#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse

parser = argparse.ArgumentParser(description="Remove charge column from LAMMPS data file")
parser.add_argument("--input_file", "-i", type=str, default="MgO.lmp",  help="input filename")
parser.add_argument("--output_file", "-o", type=str, default="MgO_atomic.lmp",  help="output filename")
args = parser.parse_args()

def modify_lammps_file(file_path, output_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(output_path, 'w') as new_file:
        atoms_section = False
        for line in lines:
            if 'Atoms # charge' in line:
                new_file.write('Atoms # atomic\n')
                atoms_section = True
                continue

            if atoms_section and line.strip() and not line.strip().startswith('#'):
                parts = line.split()
                if len(parts) > 3:
                    new_line = ' '.join(parts[:2] + parts[3:]) + '\n'
                    new_file.write(new_line)
                else:
                    new_file.write(line)
            else:
                new_file.write(line)

modify_lammps_file(args.input_file, args.output_file)
