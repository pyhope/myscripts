#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dpdata import LabeledSystem
import numpy as np
import glob
import argparse

parser = argparse.ArgumentParser(description="Calculate the RMSE between VASP and NN forces")
parser.add_argument("--input_directory", "-i", default = '.', type=str,  help="Input directory that contains OUTCAR and LAMMPS dump files")

args = parser.parse_args()

def extract_outcar(outcar):
    ls     = LabeledSystem(outcar,fmt='outcar') 
    forces = ls['forces']
    return forces[0]

def extract_dump(dump):
    forces = []
    with open(dump, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "ITEM: ATOMS" in line:
                # start of forces data
                for j in range(i+1, len(lines)):
                    data = lines[j].split()
                    fx, fy, fz = float(data[8]), float(data[9]), float(data[10])
                    forces.append((fx, fy, fz))
    return np.array(forces)

def RMSE(x1,x2):
    dx = x1-x2
    rmse = np.sqrt(np.sum(dx**2)/len(dx))
    return rmse

def dev_vasp_nn(vasp, nn, natoms=161):
    vasp_force = vasp.reshape(len(vasp)*3) # each step, there are natoms*3 forces
    nn_force = nn.reshape(len(nn)*3)
    rmse_f = RMSE(vasp_force, nn_force)

    return rmse_f

path = args.input_directory
forces_vasp = extract_outcar(path + '/OUTCAR')
forces_nn = extract_dump(glob.glob(path + '/*.dump')[0])
rmse = dev_vasp_nn(forces_vasp, forces_nn)
print('RMSE = ', rmse)
with open(path + '/RMSE_F', 'w') as f:
    f.write(str(rmse))

# np.savetxt('forces_vasp.txt', forces_vasp)
# np.savetxt('forces_nn.txt', forces_nn)
