#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

parser = argparse.ArgumentParser(description="cut the lammps dump file")
parser.add_argument("--input_file", "-i", type=str, default="nvt.dump",  help="original dump filename")
parser.add_argument("--output_file", "-o", type=str, default="new_nvt.dump",  help="new dump filename")
parser.add_argument("--endstep", "-es", type=int, help="The cut timestep of the dump file (this step is not included in the new dump file)")
args = parser.parse_args()

endstep = str(args.endstep)

with open(args.input_file, 'r') as datain:
    with open(args.output_file, 'w') as dataout:
        line = datain.readline()
        dataout.write(line)
        while line:
            line = datain.readline()
            if line.rstrip() == "ITEM: TIMESTEP":
                nextline = datain.readline()
                if nextline.rstrip() == endstep:
                    break
                else:
                    dataout.write(line)
                    dataout.write(nextline)
            else:
                dataout.write(line)
