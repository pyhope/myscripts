#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import my_pyplot as mpt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str,  help="Input file")
args = parser.parse_args()

dirctions = ['x', 'y', 'z']
c = ['C0', 'C1', 'C2']

with open(args.input_file, 'rb') as file:
    T, P, tcorr, jxyz, kxyz = pkl.load(file)
tcorr_idx = range(len(tcorr))

fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)

ax[0].plot(tcorr, np.mean(jxyz, axis=0))
ax[1].plot(tcorr, np.mean(kxyz, axis=0), c='k')
# for i, k in enumerate(kxyz):
#     ax[1].plot(tcorr, k, c=c[i], label=dirctions[i])

ax[1].set_xlabel('t (ps)')
ax[1].legend(fancybox=False, edgecolor='black')
#plt.gca().set_xlim(2e-4, 100)

for a in ax:
    a.set_xscale('log')
    a.set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')

plt.show()