#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys

parser = argparse.ArgumentParser(description="Extract frames from LAMMPS dump file.")
parser.add_argument("--input_filename", "-i", type=str, default='nvt.dump', help="Input filename (default: nvt.dump)")
parser.add_argument("--frequency", "-f", type=int, default=1000000, help="Frequency of frames to be extracted (default: 1000000)")
parser.add_argument("--start", "-s", type=int, default=0, help="Starting timestep number (default: 0)")
parser.add_argument("--start_from_last", "-sl", type=int, help="Start timestep number (count from the last frame, default: None)")
parser.add_argument("--logfilename", "-l", type=str, default='log.lammps', help="LAMMPS log filename (default: log.lammps)")
parser.add_argument("--combined", "-c", action='store_true', help="Combine all extracted frames into a single output file (default: False)")
parser.add_argument("--combined_filename", type=str, default='extracted.dump', help="Output filename if using --combined (default: extracted.dump)")
parser.add_argument("--sequential", "-q", action='store_true', help="Ignore TIMESTEP number and assign frame indices sequentially starting from 1")

args = parser.parse_args()

# Sanity check
if args.sequential and args.start_from_last:
    sys.exit("Error: --sequential (-q) mode is incompatible with --start_from_last (-sl).")

# Get true start timestep if not using sequential mode
if args.start_from_last and not args.sequential:
    from lammps_logfile import File
    log = File(args.logfilename)
    steps = log.get("Step")
    start = steps[-1] - args.start_from_last
else:
    start = args.start

datain = open(args.input_filename, 'r')

if args.combined:
    combined_out = open(args.combined_filename, 'w')

# Read frames
line = datain.readline()
frame_index = 0  # for sequential mode

while line:
    if "ITEM: TIMESTEP" in line:
        timestep_line = datain.readline()
        original_timestep = int(timestep_line.strip())
        frame_index += 1

        # Determine whether to output this frame
        if args.sequential:
            should_output = (frame_index % args.frequency == 0)
            output_id = frame_index // args.frequency
        else:
            should_output = (original_timestep >= start and original_timestep % args.frequency == 0)
            output_id = original_timestep // args.frequency

        if should_output:
            if args.combined:
                dataout = combined_out
            else:
                output_filename = f"{output_id}.dump"
                dataout = open(output_filename, 'w')

            # Overwrite timestep number if in sequential mode
            dataout.write("ITEM: TIMESTEP\n")
            if args.sequential:
                dataout.write(f"{frame_index}\n")
            else:
                dataout.write(timestep_line)

            # Write rest of the frame
            while True:
                line = datain.readline()
                if not line or "ITEM: TIMESTEP" in line:
                    break
                dataout.write(line)

            if not args.combined:
                dataout.close()
            else:
                continue  # avoid extra readline below

    line = datain.readline()

datain.close()
if args.combined:
    combined_out.close()
