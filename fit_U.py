#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy.optimize import curve_fit

def linear_func(x, m, b):
    return m * x + b

parser = argparse.ArgumentParser(description="Fit Hubbard U parameter")
parser.add_argument("--input_file", "-i", type=str, default="n_elec.txt",  help="input file name")
parser.add_argument("--output_file", "-o", type=str, default="U.txt",  help="output file name")

args = parser.parse_args()

V, n_nscf, n_scf = np.loadtxt(args.input_file, unpack=True)

params, cov = curve_fit(linear_func, V, n_nscf)
slope_nscf, _ = params

params, cov = curve_fit(linear_func, V, n_scf)
slope_scf, _ = params

U = 1/slope_scf - 1/slope_nscf

with open(args.output_file, "w") as f:
    f.write(f"{U:.2f}\n")

print(f"U = {U:.2f}")
