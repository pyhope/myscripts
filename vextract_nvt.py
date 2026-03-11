#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

parser = argparse.ArgumentParser(description="Extract configurations from an XDATCAR file")
parser.add_argument("--input_file", "-i", type=str, default="./XDATCAR", help="Path to the XDATCAR file")
parser.add_argument("--output_dir", "-o", type=str, default="./frames", help="Directory for the output POSCAR files")
parser.add_argument("--frequency", "-f", type=int, default=100, help="Frequency of configurations to extract")
parser.add_argument("--num_frames", "-n", type=int, default=5, help="Number of configurations to extract")
parser.add_argument(
    "--start_frame", "-s",
    type=int,
    default=None,
    help="Start frame number for initial mode (1-based indexing). If not provided, extraction is from the end."
)
args = parser.parse_args()


def extract_frames_from_xdatcar(xdatcar_file, output_dir, interval, num_frames, start_frame=None):
    if interval <= 0:
        raise ValueError("frequency/interval must be a positive integer.")
    if num_frames <= 0:
        raise ValueError("num_frames must be a positive integer.")

    with open(xdatcar_file, "r") as file:
        lines = file.readlines()

    header = lines[:7]
    num_atoms_line = lines[6]

    config_indices = [i for i, line in enumerate(lines) if line.startswith("Direct configuration=")]
    total_configs = len(config_indices)

    if total_configs == 0:
        raise ValueError("No 'Direct configuration=' entries found in the XDATCAR file.")

    if start_frame is not None:
        if start_frame < 1:
            raise ValueError("start_frame must be >= 1.")
        if start_frame > total_configs:
            print(
                f"Warning: start_frame ({start_frame}) is larger than total number of frames "
                f"({total_configs}). No frame will be extracted."
            )
            return 0

    selected_configs = []

    if start_frame is None:
        # Extract from the end
        for i in range(total_configs - 1, -1, -interval):
            selected_configs.append(config_indices[i])
            if len(selected_configs) == num_frames:
                break
    else:
        # Extract from given start frame (1-based indexing)
        current_frame = start_frame - 1
        while current_frame < total_configs and len(selected_configs) < num_frames:
            selected_configs.append(config_indices[current_frame])
            current_frame += interval

    actual_num_frames = len(selected_configs)

    if actual_num_frames < num_frames:
        print(
            f"Warning: requested {num_frames} frames with interval {interval}, "
            f"but only {actual_num_frames} frame(s) are available from the given start point."
        )

    num_atoms = sum(map(int, num_atoms_line.split()))

    for idx, config_index in enumerate(selected_configs):
        frame_lines = header + ["Direct\n"]
        frame_lines += lines[config_index + 1: config_index + 1 + num_atoms]

        config_num = lines[config_index].split()[-1]
        poscar_filename = os.path.join(output_dir, f"POSCAR.{idx + 1}")

        with open(poscar_filename, "w") as poscar_file:
            poscar_file.writelines(frame_lines)

        print(f"Configuration {config_num} written to {poscar_filename}")

    print(f"Requested {num_frames} frame(s); actually extracted {actual_num_frames} frame(s).")
    return actual_num_frames


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

extract_frames_from_xdatcar(
    args.input_file,
    args.output_dir,
    args.frequency,
    args.num_frames,
    start_frame=args.start_frame
)