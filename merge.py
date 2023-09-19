#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description="Merge the dump files from different runs. ENDSTEPS is required.")
parser.add_argument("dir_list", nargs="+", type=str, help="list of directories")
parser.add_argument("--input_filename", "-i", type=str, default="nvt.dump",  help="input dump filename")
parser.add_argument("--output_filename", "-o", type=str, default="full_nvt.dump",  help="output dump filename")

args = parser.parse_args()
dirnames = args.dir_list
f_in = args.input_filename
f_out = args.output_filename

subprocess.run(["findend"] + dirnames)
print("ENDSTEPS file generated!")

endsteps = np.loadtxt('ENDSTEPS', unpack=True, dtype=str)
dataout = open(f_out, 'w')
for i in range(len(dirnames)):
    datain = open(dirnames[i]+'/' + f_in, 'r')
    print("Processing " + dirnames[i] + '/' + f_in + " ...")
    line = datain.readline()
    dataout.write(line)
    while line:
        line = datain.readline()
        if line.rstrip() == "ITEM: TIMESTEP" and i < len(dirnames)-1:
            nextline = datain.readline()
            if nextline.rstrip() == endsteps[i]:
                break
            else:
                dataout.write(line)
                dataout.write(nextline)
        else:
            dataout.write(line)
    datain.close()
dataout.close()