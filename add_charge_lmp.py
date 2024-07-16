#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description='Add charge to lammps data file')
parser.add_argument("--input_file", "-i", type=str, default="MgFeO.lmp",  help="original data filename")
parser.add_argument("--output_file", "-o", type=str, default="charge.lmp",  help="new data filename")
parser.add_argument("--isMgFeO", "-Fe", default=True, action='store_false', help="Defualt: this is MgFeO")
args = parser.parse_args()

def assign_charge(lammps_data_file, output_file):
    with open(lammps_data_file, 'r') as file:
        lines = file.readlines()

    output_data = []
    for line in lines:
        if "Atoms # atomic" in line:
            output_data.append("Atoms # charge\n\n")
            break
        else:
            output_data.append(line)
    
    atoms_lines = lines[len(output_data):]
    for atom in atoms_lines:
        tokens = atom.strip().split()
        if len(tokens) <= 1:
            continue
        atom_id = int(tokens[0])
        atom_type = int(tokens[1])
        x = tokens[2]
        y = tokens[3]
        z = tokens[4]
        
        if args.isMgFeO:
            if atom_type == 1:
                charge = 1.2
            elif atom_type == 2:
                charge = 1.2
            elif atom_type == 3:
                charge = -1.2
            else:
                charge = 0
        else:
            if atom_type == 1:
                charge = 1.2
            elif atom_type == 2:
                charge = -1.2
            else:
                charge = 0
        
        modified_atom_line = f"{atom_id} {atom_type} {charge} {x} {y} {z}\n"
        output_data.append(modified_atom_line)

    with open(output_file, 'w') as file:
        file.writelines(output_data)

assign_charge(args.input_file, args.output_file)
