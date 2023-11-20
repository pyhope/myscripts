#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description="Extract frames from LAMMPS dump file.")
parser.add_argument("--input_filename", "-i", type=str, default='nvt.dump', help="Input filename (default: nvt.dump)")

args = parser.parse_args()

datain = open(args.input_filename, 'r')

# Find and record the target frames
line = datain.readline()

while line:
    if "ITEM: TIMESTEP" in line:
        timestep_line = datain.readline()
        thisframe = int(timestep_line.strip())
        
        # Check if current frame is divisible by 1000000
        if thisframe % 1000000 == 0:
            output_filename = str(thisframe // 1000000) + '.dump'
            with open(output_filename, 'w') as dataout:
                dataout.write("ITEM: TIMESTEP\n")
                dataout.write(timestep_line)
                while True:  # Write the rest of the frame data until next "ITEM: TIMESTEP" or EOF
                    line = datain.readline()
                    if not line or "ITEM: TIMESTEP" in line:
                        break
                    dataout.write(line)
    
    line = datain.readline()

datain.close()
