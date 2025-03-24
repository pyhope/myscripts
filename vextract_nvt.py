#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

parser = argparse.ArgumentParser(description="Extract configurations from an XDATCAR file")
parser.add_argument("--input_file", "-i", type=str, default="./XDATCAR",  help="Path to the XDATCAR file")
parser.add_argument("--output_dir", "-o", type=str, default="./frames", help="Directory for the output POSCAR files")
parser.add_argument("--frequency", "-f", type=int, default=100, help="Frequency of configurations to extract")
parser.add_argument("--num_frames", "-n", type=int, default=5, help="Number of configurations to extract")
parser.add_argument("--start_frame", "-s", type=int, default=None, help="Start frame number for initial mode (1-based indexing). If not provided, extraction is from the end.")
args = parser.parse_args()


def extract_frames_from_xdatcar(xdatcar_file, output_dir, interval, num_frames, start_frame=None):
    with open(xdatcar_file, 'r') as file:
        lines = file.readlines()

    header = lines[:7]
    num_atoms_line = lines[6]

    config_indices = [i for i, line in enumerate(lines) if line.startswith("Direct configuration=")]

    total_configs = len(config_indices)
    if total_configs == 0:
        raise ValueError("No 'Direct configuration' found.")

    selected_configs = []

    if start_frame is None:
        for i in range(total_configs - 1, -1, -interval):
            selected_configs.append(config_indices[i])
            if len(selected_configs) == num_frames:
                break
    else:
        current_frame = start_frame - 1
        while current_frame < total_configs and len(selected_configs) < num_frames:
            selected_configs.append(config_indices[current_frame])
            current_frame += interval

    if len(selected_configs) < num_frames:
        raise ValueError(f"Insufficient configurations to extract {num_frames} frames with interval {interval} from given start point.")

    num_atoms = sum(map(int, num_atoms_line.split()))
    for idx, config_index in enumerate(selected_configs):
        frame_lines = header + ['Direct\n']
        frame_lines += lines[config_index + 1: config_index + 1 + num_atoms]

        config_num = lines[config_index].split()[-1]
        poscar_filename = f"{output_dir}/POSCAR.{idx + 1}"
        with open(poscar_filename, 'w') as poscar_file:
            poscar_file.writelines(frame_lines)
        print(f"Configuration {config_num} written to {poscar_filename}")


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

extract_frames_from_xdatcar(
    args.input_file,
    args.output_dir,
    args.frequency,
    args.num_frames,
    start_frame=args.start_frame
)
