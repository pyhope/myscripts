#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ase.io import read, write
import argparse

parser = argparse.ArgumentParser(description="Convert LAMMPS dump file to POSCAR with velocities")
parser.add_argument("--input_file", "-i", type=str, default="selected.dump",  help="input .dump file")

args = parser.parse_args()
dump_file = args.input_file

frame = read(dump_file, index=0, format="lammps-dump-text")

chemical_symbols = []

for atom in frame:
    if atom.symbol == 'H':
        chemical_symbols.append('Mg')
    elif atom.symbol == 'He':
        chemical_symbols.append('Si')
    elif atom.symbol == 'Li':
        chemical_symbols.append('O')
    elif atom.symbol == 'Be':
        chemical_symbols.append('H')

frame.set_chemical_symbols(chemical_symbols)

write('POSCAR', frame, format='vasp')

def read_lammps_dump(filename):
    velocities = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "ITEM: ATOMS" in line:
                # start of velocities data
                for j in range(i+1, len(lines)):
                    data = lines[j].split()
                    vx, vy, vz = float(data[5])/1000, float(data[6])/1000, float(data[7])/1000  # Divide by 1000 to convert to angstroms/femtosecond
                    velocities.append((vx, vy, vz))
    return velocities

def append_to_poscar(velocities, poscar_file="POSCAR"):
    with open(poscar_file, 'a') as f:
        f.write("\n")
        for v in velocities:
            f.write("{:.16f} {:.16f} {:.16f}\n".format(v[0], v[1], v[2]))

velocities = read_lammps_dump(dump_file)
append_to_poscar(velocities)