#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Trim the initial frames of an Extended XYZ file")
parser.add_argument("--input_filename", "-i", type=str, default='ASAP-desc.xyz', help="Input filename (default: ASAP-desc.xyz)")
parser.add_argument("--output_filename", "-o", type=str, default='ASAP-desc-trimmed.xyz', help="Output filename (default: ASAP-desc-trimmed.xyz)")
parser.add_argument("--num_frames_to_remove", "-n", type=int, default=500, help="Number of initial frames to remove (default: 500)")

args = parser.parse_args()
n = args.num_frames_to_remove

with open(args.input_filename, 'r') as infile, open(args.output_filename, 'w') as outfile:
    frame_count = 0
    while frame_count < n:
        atom_count_line = infile.readline()
        if not atom_count_line:
            break
        atom_count = int(atom_count_line.strip())
        
        infile.readline()
        for _ in range(atom_count):
            infile.readline()
        frame_count += 1

    while True:
        line = infile.readline()
        if not line:
            break
        outfile.write(line)
