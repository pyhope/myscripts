#!/usr/bin/env python3

import numpy as np
import lammps_logfile
import argparse

parser = argparse.ArgumentParser(description="calculate properties from lammps log files")
parser.add_argument("-i", "--input_file", type=str, default="log.lammps",  help="LAMMPS output file")
parser.add_argument("-p", "--properties",  nargs="+", type=str, help="Properties to calculate")
parser.add_argument("-o", "--out_file",type=str,default='out.txt', help="Out file name")
parser.add_argument("-ss", "--start_step",type=int,default=5000, help="Consider the properties from this step")
args = parser.parse_args()

log = lammps_logfile.File(args.input_file)
f = args.out_file
s = args.start_step

if args.properties:
    props = args.properties
else:
    props = ["Temp", "Press", "Lx", "Ly", "Lz", "Volume"]

Data = dict()
for i in props:
    Data[i] = dict()
    Data[i]['avg'] = np.mean(log.get(i)[s:])
    Data[i]['std'] = np.std(log.get(i)[s:])

with open(args.out_file, 'w') as file:
    file.write('#' + ' '.join(props) + '\n')
    file.write(' '.join([str(Data[i]['avg']) for i in props]) + '\n')
    file.write(' '.join([str(Data[i]['std']) for i in props]) + '\n')