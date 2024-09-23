#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Calculate volume of POSCAR")
parser.add_argument("--input_file", "-i", type=str, default="POSCAR",  help="input file name")
parser.add_argument("--output_file", "-o", type=str, default="vol.txt",  help="output file name")

args = parser.parse_args()

import numpy as np

def calculate_parallelepiped_volume(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    vector_a = np.array([float(x) for x in lines[2].split()])
    vector_b = np.array([float(x) for x in lines[3].split()])
    vector_c = np.array([float(x) for x in lines[4].split()])

    cross_product = np.cross(vector_b, vector_c)
    volume = abs(np.dot(vector_a, cross_product))

    return volume

file_path = args.input_file
volume = calculate_parallelepiped_volume(file_path)
print(f"Volume: {volume:.2f}")
with open(args.output_file, 'w') as file:
    file.write(f"{volume:.2f}\n")
