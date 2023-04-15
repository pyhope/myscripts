#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Oct  9 12:25:43 2022
@author: Yihang Peng
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_filename","-i",type=str, default = 'npt.dump' ,help="input filename (npt.dump)")
parser.add_argument("--output_filename","-o",type=str, default = 'selected.dump' ,help="output filename (selected.dump)")
parser.add_argument("--timestep","-t",type=str, default = '50000' ,help="target frame")

args   = parser.parse_args()

timestep = args.timestep

datain = open(args.input_filename, 'r')

# Find and record the target frame
line = datain.readline()
while line:
    while line.rstrip() != "ITEM: TIMESTEP":
        if line:
            line = datain.readline()
        else:
            print("Not Found!")
            break
    else:
        writedata = line
        line = datain.readline()
        if line.rstrip() == timestep:
            writedata += line
            break
writedata += datain.readline()
line = datain.readline()
atom_num = int(line)
writedata += line
while line.split(' ')[0] != str(atom_num):
    if line:
        line = datain.readline()
        writedata += line
    else:
        print("Error!")

# Write data to new file
dataout = open(args.output_filename, 'w')
dataout.write(writedata)