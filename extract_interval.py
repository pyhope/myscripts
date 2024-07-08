#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description="Extract frames from LAMMPS dump file.")
parser.add_argument("--input_filename", "-i", type=str, default='nvt.dump', help="Input filename (default: nvt.dump)")
parser.add_argument("--frequency", "-f", type=int, default=1000000, help="Frequency of frames to be extracted (default: 1000000)")
parser.add_argument("--start", "-s", type=int, default=0, help="Starting timestep number (default: 0)")
parser.add_argument("--start_from_last", "-sl", type=int, help="Start timestep number (count from the last frame, default: None)")
parser.add_argument("--logfilename", "-l", type=str, default='log.lammps', help="LAMMPS log filename (default: log.lammps)")

args = parser.parse_args()

if args.start_from_last:
    from lammps_logfile import File
    log = File(args.logfilename)
    steps = log.get("Step")
    start = steps[-1] - args.start_from_last
else:
    start = args.start

datain = open(args.input_filename, 'r')

# Find and record the target frames
line = datain.readline()

while line:
    if "ITEM: TIMESTEP" in line:
        timestep_line = datain.readline()
        thisframe = int(timestep_line.strip())
        
        # Check if current frame is divisible by args.frequency
        if thisframe >= start and thisframe % args.frequency == 0:
            output_filename = str(thisframe // args.frequency) + '.dump'
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
