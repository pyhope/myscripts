#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Yihang Peng

import MDAnalysis as mda
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Calculate the average distance between helium atoms")
parser.add_argument("--input_file", "-i", type=str, default="tmp.pdb",  help="pdb file of the crystal structure")
args = parser.parse_args()

u = mda.Universe(args.input_file)

atoms = u.select_atoms("name HE")

with open("INDEX", "w") as f:
    for atom in atoms:
        f.write(str(atom.index + 1) + "\n")

distances = []
for i in range(len(atoms)):
    i_distances = []
    for j in range(len(atoms)):
        if i == j:
            continue
        distance = mda.lib.distances.distance_array(atoms[i].position[None, :], atoms[j].position[None, :], box=u.dimensions)[0, 0]
        i_distances.append(distance)
    distances.append(np.min(i_distances))

average_distance = np.mean(distances)

with open("DISTANCE", "w") as f:
    f.write(str(average_distance) + "\n")