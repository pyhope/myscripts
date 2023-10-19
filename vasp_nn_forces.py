#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dpdata import LabeledSystem
import numpy as np
import glob

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
                    fx, fy, fz = float(data[8]), float(data[9]), float(data[10])  # Divide by 1000 to convert to angstroms/femtosecond
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

path = '/scratch/gpfs/yp0007/project/diffusion/t2.5p25-81/en/582810/'
forces_vasp = extract_outcar(path + 'OUTCAR')
forces_nn = extract_dump(glob.glob(path + '*.dump')[0])
rmse = dev_vasp_nn(forces_vasp, forces_nn)
print('rmse = ', rmse)

# np.savetxt('forces_vasp.txt', forces_vasp)
# np.savetxt('forces_nn.txt', forces_nn)
