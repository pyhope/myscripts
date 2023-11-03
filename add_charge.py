#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def assign_charge(lammps_dump_file, output_file):
    with open(lammps_dump_file, 'r') as file:
        lines = file.readlines()

    output_data = []
    for line in lines:
        if "ITEM: ATOMS id type x y z" in line:
            output_data.append("ITEM: ATOMS id type q x y z\n")
        else:
            output_data.append(line)
        
        if "ITEM: ATOMS" in line:
            break
    
    atoms_lines = lines[len(output_data):]
    for atom in atoms_lines:
        tokens = atom.strip().split()
        atom_id = int(tokens[0])
        atom_type = int(tokens[1])
        x = tokens[2]
        y = tokens[3]
        z = tokens[4]
        
        if atom_type == 1:
            charge = 1.2
        elif atom_type == 2:
            charge = -1.2
        else:
            charge = 0  # or any default value you want for other types
        
        modified_atom_line = f"{atom_id} {atom_type} {charge} {x} {y} {z}\n"
        output_data.append(modified_atom_line)

    with open(output_file, 'w') as file:
        file.writelines(output_data)

lammps_dump_file = 'test.dump'
output_file = 'test_charge.dump'
assign_charge(lammps_dump_file, output_file)
