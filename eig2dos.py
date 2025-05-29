#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--smear", "-s", type=float, default=0.35355, help="Smearing width for DOS calculation")
args = parser.parse_args()

smear = args.smear

# --- Read Efermi from DOSCAR ---
with open('DOSCAR', 'r') as f:
    for _ in range(5):
        f.readline()
    linfo = f.readline()
    Emax, Emin, Nbins, Efermi, _ = np.asarray(linfo.split(), dtype=float)

# --- Determine ISPIN from INCAR ---
ISPIN = 1
if os.path.exists('INCAR'):
    with open('INCAR', 'r') as f:
        for line in f:
            if 'ISPIN' in line:
                ISPIN = int(line.split('=')[1].strip())
                break

# --- Read EIGENVAL ---
col_names = ('kptind', 'kx', 'ky', 'kz', 'kw', 'energy', 'occupancy')
with open('EIGENVAL', 'r') as f:
    for _ in range(5):
        f.readline()
    Nele, Nkpts, Nbands = np.asarray(f.readline().split(), dtype=int)
    f.readline()

    kpts = np.zeros((Nkpts, 4))

    if ISPIN == 2:
        keconcat_up = np.zeros((Nkpts, Nbands, 7))
        keconcat_down = np.zeros((Nkpts, Nbands, 7))
    else:
        keconcat = np.zeros((Nkpts, Nbands, 7))

    for i in range(Nkpts):
        kpts[i] = np.asarray(f.readline().split(), dtype=float)
        kptstl = np.tile(kpts[i], (Nbands, 1))

        if ISPIN == 2:
            eigens_up = np.zeros((Nbands, 2))
            eigens_down = np.zeros((Nbands, 2))

            for n in range(Nbands):
                parts = f.readline().split()
                if len(parts) < 5:
                    raise ValueError(f"Invalid EIGENVAL line format at band {n}, k-point {i}")
                eigens_up[n] = [float(parts[1]), float(parts[3])]
                eigens_down[n] = [float(parts[2]), float(parts[4])]

            keconcat_up[i] = np.concatenate((np.ones((Nbands, 1)) * (i + 1), kptstl, eigens_up), axis=1)
            keconcat_down[i] = np.concatenate((np.ones((Nbands, 1)) * (i + 1), kptstl, eigens_down), axis=1)
        else:
            eigens = np.zeros((Nbands, 2))
            for n in range(Nbands):
                parts = f.readline().split()
                eigens[n] = np.asarray(parts[1:3], dtype=float)
            keconcat[i] = np.concatenate((np.ones((Nbands, 1)) * (i + 1), kptstl, eigens), axis=1)

        f.readline()  # blank line after each k-point

# --- Store results ---
if ISPIN == 2:
    kedf_up = pd.DataFrame(np.reshape(keconcat_up, (Nbands * Nkpts, 7)), columns=col_names)
    kedf_down = pd.DataFrame(np.reshape(keconcat_down, (Nbands * Nkpts, 7)), columns=col_names)
    kedf_up.to_csv('ke_summary_up.csv', index=False)
    kedf_down.to_csv('ke_summary_down.csv', index=False)
else:
    kedf = pd.DataFrame(np.reshape(keconcat, (Nbands * Nkpts, 7)), columns=col_names)
    kedf.to_csv('ke_summary.csv', index=False)

# --- DOS calculation helpers ---
def gaussian(x, mu, sig, I):
    return np.exp(-(x - mu)**2. / (2. * sig**2.)) * I / (sig * np.sqrt(2. * np.pi))

def y1new(x, x1, y1, smear):
    return np.sum([gaussian(x, x1[i], smear, y1[i]) for i in range(x1.size)], axis=0)

def compute_dos(kedf, tag=''):
    emin = kedf['energy'].min()
    emax = kedf['energy'].max()
    ebin = np.arange(emin, emax, 0.1)
    dos = np.zeros(ebin.size - 1)
    dos_occ = np.zeros(ebin.size - 1)

    for i in range(1, Nkpts + 1):
        kei = kedf[kedf['kptind'] == i]
        e = kei['energy']
        occup = kei['occupancy']
        kw = kei.iloc[0]['kw']
        dos += kw * np.histogram(e, bins=ebin)[0]
        dos_occ += kw * np.histogram(e, bins=ebin, weights=occup)[0]

    ebin_sm = np.arange(emin - 10 * smear, emax + 10 * smear, smear)
    dossm = y1new(ebin_sm, ebin[:-1], dos, smear)
    dossm_occ = y1new(ebin_sm, ebin[:-1], dos_occ, smear)

    np.savetxt(f'dos{tag}.txt', np.column_stack((ebin[:-1] - Efermi, dos)), fmt='%f', delimiter='\t ',
               header=f'# E-Ef(eV)\t dos{tag}', comments='')
    np.savetxt(f'dos_occup{tag}.txt', np.column_stack((ebin[:-1] - Efermi, dos_occ)), fmt='%f', delimiter='\t ',
               header=f'# E-Ef(eV)\t dos_occ{tag}', comments='')
    np.savetxt(f'dos_sm{tag}.txt', np.column_stack((ebin_sm - Efermi, dossm)), fmt='%f', delimiter='\t ',
               header=f'# E-Ef(eV)\t dossm{tag}', comments='')
    np.savetxt(f'dos_occup_sm{tag}.txt', np.column_stack((ebin_sm - Efermi, dossm_occ)), fmt='%f', delimiter='\t ',
               header=f'# E-Ef(eV)\t dossm_occ{tag}', comments='')

    return ebin[:-1] - Efermi, dos, dos_occ, ebin_sm - Efermi, dossm, dossm_occ

# --- Run and plot ---
if ISPIN == 2:
    x1, dos1, dos1_occ, x1sm, dossm1, dossm1_occ = compute_dos(kedf_up, '_up')
    x2, dos2, dos2_occ, x2sm, dossm2, dossm2_occ = compute_dos(kedf_down, '_down')
else:
    x, dos, dos_occ, xsm, dossm, dossm_occ = compute_dos(kedf)
