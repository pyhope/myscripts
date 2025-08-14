#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("-r", "--run_dir", type=str, default="../1_run", help="Phq running directory")
parser.add_argument("-gv", "--gv_dir", type=str, default="../3_gv", help="Group velocity calculation directory")
parser.add_argument("-c", "--frequency_cutoff", type=float, default=400, help="Cutoff frequency (cm^-1)")
args = parser.parse_args()

THz_to_cm_inv = 33.35641
cut = args.frequency_cutoff

def power_law(x, a, b):
    return a * x**b

def fit_and_plot(x, y, color, linestyle):
    mask = ~np.isnan(x) & ~np.isnan(y)
    x, y = x[mask], y[mask]

    popt, _ = curve_fit(power_law, x, y, p0=[1e-6, 5])
    a, b = popt
    xfit = np.linspace(x.min(), x.max(), 200)
    yfit = power_law(xfit, a, b)
    axes[0].plot(xfit, yfit, linestyle, color=color, lw=1)
    axes[1].plot(xfit, 1 / (yfit / THz_to_cm_inv * 2), linestyle, color=color, lw=1)
    print(f"Fitted parameters: a={a:.2e}, b={b:.2f}")

    data = np.loadtxt(f"{args.gv_dir}/group_velocity.txt")
    freq_fit = data[:,1]
    tau_fit = 1 / (power_law(freq_fit, a, b) / THz_to_cm_inv * 2)
    tau_fit[0:3] = 0
    np.savetxt(f"tau_fit.dat", np.column_stack((freq_fit, tau_fit)), fmt='%.6f', header='Frequency(cm^-1) Tau(ps)')

fig, axes = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(12, 6))
fig.subplots_adjust(wspace=0.2)

data = np.loadtxt(f"{args.run_dir}/tau_fit.tau")
linewidth = data[:,1][3:]
tau_ps = 1 / (linewidth / THz_to_cm_inv * 2)
tau_ps = np.where(tau_ps > 5, np.nan, tau_ps)
data = np.loadtxt(f"{args.run_dir}/frequency.freq")
freq_cm = data[:,2][3:]

axes[0].scatter(freq_cm, linewidth, s=8, ec='C0', fc='none')
axes[1].scatter(freq_cm, tau_ps, s=8, ec='C0', fc='none')

mask_low = freq_cm < cut
mask_high = freq_cm >= cut

omega_low = freq_cm[mask_low]
invtau_low = 1 / tau_ps[mask_low]
omega_high = freq_cm[mask_high]
invtau_high = 1 / tau_ps[mask_high]

if len(omega_low) > 2:
    fit_and_plot(omega_low, invtau_low * THz_to_cm_inv / 2, 'C0', '--')
if len(omega_high) > 2:
    fit_and_plot(omega_high, invtau_high * THz_to_cm_inv / 2, 'C0', '-.')

axes[0].set_ylabel('Phonon linewidth (cm$^{-1}$)')
axes[1].set_ylabel('Phonon lifetime (ps)')
axes[1].set_ylim(0, 3)

for ax in axes:
    ax.minorticks_on()
    ax.set_xlabel('Frequency (cm$^{-1}$)')

plt.savefig(f"fit_tau.png", dpi=300, bbox_inches='tight')
