#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_filename", type=str, help="input filename")
parser.add_argument("--output_filename","-o",type=str, default = 'test.txt' , help="output filename")

args = parser.parse_args()

def extract_values(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    results = []

    for i, line in enumerate(lines):
        if 'number of test data' in line:
            for j in range(2, 5):
                target_line = lines[i + j].strip()
                values = target_line.split()
                if len(values) > 1:
                    results.append(values[-2])

    with open(output_file, 'w') as f:
        for value in results:
            f.write(value + '\n')
        f.write('\n')

input_file = args.input_filename
output_file = args.output_filename

extract_values(input_file, output_file)
