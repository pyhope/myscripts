#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse

# Argument parser for command line options
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--energy', help='Energy file', required=True)
parser.add_argument('-f', '--force', help='Force file', required=True)
parser.add_argument('-v', '--virial', help='Virial file', required=True)
args = parser.parse_args()

def rmse(org,pred): # same as dp_test l2err
    dif = pred - org
    return np.sqrt(np.mean(dif**2))

e=np.loadtxt(args.energy)
f=np.loadtxt(args.force)
v=np.loadtxt(args.virial)

print("N\t", len(e))
print("e\t", rmse(e[:,0],e[:,1]))
print("f\t", rmse(f[:,:3],f[:,3:]))
print("v\t", rmse(v[:,:9],v[:,9:]))
