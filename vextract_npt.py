#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

parser = argparse.ArgumentParser(description="Extracts configurations from a VASP XDATCAR file.")
parser.add_argument("--input_file", "-i", type=str, default="./XDATCAR",  help="Path to the XDATCAR file")
parser.add_argument("--start_block", "-s", type=int, default=1000, help="Block number to start extracting configurations from")
parser.add_argument("--frequency", "-f", type=int, default=100, help="Frequency of configurations to extract")
args = parser.parse_args()

def split_and_store_blocks(filename, start_block_num, frequency):
    with open(filename, 'r') as file:
        lines = file.readlines()

    block_count = 0
    current_block = []
    for line in lines:
        if line.startswith("unknown"):
            if block_count >= start_block_num and (block_count - start_block_num) % frequency == 0:
                save_block(current_block, block_count)
            current_block = [line]
            block_count += 1
        else:
            current_block.append(line)

    # Save the last block if it meets the criteria
    if block_count >= start_block_num and (block_count - start_block_num) % frequency == 0:
        save_block(current_block, block_count)

def save_block(block, block_num):
    block_dir = f"{block_num}"
    os.makedirs(block_dir, exist_ok=True)
    with open(os.path.join(block_dir, "POSCAR"), 'w') as file:
        file.writelines(block)

split_and_store_blocks(args.input_file, args.start_block, args.frequency)
