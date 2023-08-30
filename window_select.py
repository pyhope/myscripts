#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import argparse
from matplotlib import rcParams

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = '16'
rcParams['font.sans-serif'] = 'Arial'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str,  help="Input file")
args = parser.parse_args()

with open(args.input_file, 'rb') as file:
    T, P, tcorr, jxyz, kxyz = pkl.load(file)
tcorr_idx = range(len(tcorr))

fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)

ax[0].plot(tcorr, np.mean(jxyz, axis=0))
ax[1].plot(tcorr, np.mean(kxyz, axis=0))
ax[1].set_xlabel('t (ps)')
ax[1].set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')

plt.gca().set_xlim(2e-4, 100)

ax[1].set_xscale('log')
plt.show()