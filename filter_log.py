#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description="Convert LAMMPS dump file to POSCAR with velocities")
parser.add_argument("--input_file", "-i", type=str, default="log.lammps",  help="input log file")
parser.add_argument("--output_file", "-o", type=str, default="filtered_output.txt",  help="output file")
parser.add_argument("--beginstep","-tb",type=int, default=582800, help="target begin frame")
parser.add_argument("--endstep","-te",type=int, default=583300, help="target end frame")

args = parser.parse_args()

with open(args.input_file, 'r') as infile, open(args.output_file, 'w') as outfile:
    for line in infile:
        words = line.split()
        
        if len(words) < 2:
            continue

        first_word = words[0]
        
        try:
            value = int(first_word)
            
            if args.beginstep <= value <= args.endstep:
                outfile.write(f"{first_word} {words[1]}\n")
        except ValueError:
            continue
