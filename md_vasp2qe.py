#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser(description="Convert VASP XDATCAR to Quantum ESPRESSO md.out format")
parser.add_argument("--input_filename", "-i", type=str, default="XDATCAR",  help="input filename")
parser.add_argument("--output_filename", "-o", type=str, default="md.out",  help="output filename")
parser.add_argument("--supercell_dimension", "-d", type=int, default=2,  help="supercell dimension")
parser.add_argument("--start_step", "-s", type=int, default=2500,  help="starting step")
parser.add_argument("--end_step", "-e", type=int, default=500000,  help="ending step")

args = parser.parse_args()
xdatcar_file = args.input_filename
output_file = args.output_filename
group_count = args.supercell_dimension ** 3
start_step = args.start_step
end_step = args.end_step

def read_xdatcar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))
    natoms = sum(counts)

    # Build element list per atom
    element_list = []
    for elem, count in zip(elements, counts):
        element_list.extend([elem] * count)

    coords = []
    steps = []

    i = 7
    while i < len(lines):
        if lines[i].startswith('Direct configuration'):
            step = int(lines[i].split('=')[1])
            steps.append(step)
            i += 1
            frame_coords = []
            for _ in range(natoms):
                x, y, z = map(float, lines[i].split())
                frame_coords.append([x, y, z])
                i += 1
            coords.append(frame_coords)
        else:
            i += 1

    return elements, counts, element_list, steps, coords

def build_grouped_indices(counts, group_count):
    """
    For each element type, divide indices into equal groups.
    Return: list of index lists in the desired new order.
    """
    index_blocks = []
    start = 0
    for count in counts:
        if count % group_count != 0:
            raise ValueError(f"Atom count {count} cannot be evenly divided into {group_count} groups.")
        block_size = count // group_count
        blocks = [list(range(start + i*block_size, start + (i+1)*block_size)) for i in range(group_count)]
        index_blocks.append(blocks)
        start += count

    # Interleave across elements for each group
    reordered_indices = []
    for i in range(group_count):
        for element_blocks in index_blocks:
            reordered_indices.extend(element_blocks[i])

    return reordered_indices

def write_formatted_output(filename, element_list, steps, coords, reordered_indices, start_step, end_step):
    frame0 = coords[0]
    if end_step > len(steps):
        end_step = len(steps)
        steps = steps[start_step-1:]
        coords = coords[start_step-1:]
    else:
        steps = steps[start_step-1:end_step]
        coords = coords[start_step-1:end_step]
    with open(filename, 'w') as f:
        f.write(f"total_step = {end_step}\n\n")
        f.write("atomic_positions\n")
        for idx in reordered_indices:
            elem = element_list[idx]
            x, y, z = frame0[idx]
            f.write(f"{elem:<4}  {x:.10f}  {y:.10f}  {z:.10f}\n")
        f.write("\n")
        for step, frame in zip(steps, coords):
            f.write(f"md_step = {step}\n")
            f.write("atomic_md_positions\n")
            for idx in reordered_indices:
                elem = element_list[idx]
                x, y, z = frame[idx]
                f.write(f"{elem:<4}  {x:.10f}  {y:.10f}  {z:.10f}\n")
            f.write("\n")

elements, counts, element_list, steps, coords = read_xdatcar(xdatcar_file)
reordered_indices = build_grouped_indices(counts, group_count)
write_formatted_output(output_file, element_list, steps, coords, reordered_indices, start_step, end_step)
print(f"Conversion done. Output written to {output_file}")
