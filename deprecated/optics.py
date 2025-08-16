#!/usr/bin/python3
# -*- coding: utf-8 -*-
# author: Yihang Peng

import numpy as np
from matplotlib import pyplot as plt
import my_pyplot as mpt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", type=str, default='SIGMA', help="input file")
parser.add_argument("--output_dir", "-o", type=str, default='./', help="output directory")
parser.add_argument("--temperature", "-t", type=float, default=5000, help="Temperature in K")
parser.add_argument("--cshift", "-c", type=float, default=1e-6, help="Shift for the denominator")
args = parser.parse_args()
directory = args.output_dir
T = args.temperature

h = 6.62607015e-34 # Planck constant in J s
h_ev = 4.135667696e-15 # Planck constant in eV s
c = 299792458 # Speed of light in m/s
k_B = 1.380649e-23 # Boltzmann constant in J/K
epsilon_0 = 8.854187817e-12 # Vacuum permittivity in F/m

def kkr_single(de, eps_imag, cshift=1e-6):
    eps_imag = np.array(eps_imag)
    nedos = eps_imag.shape[0]
    cshift = complex(0, cshift)
    w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
    
    def integration_element(w_r):
        factor = w_i / (w_i**2 - w_r**2 + cshift)
        total = np.sum(eps_imag * factor)
        return total * (2/np.pi) * de + 1

    return np.real([integration_element(w_r) for w_r in w_i])

def B(nu, T):
    return 2 * h * nu**3 / c**2 / (np.exp(h * nu / (k_B * T)) - 1)

def dB_dT(nu, T):
    return (2 * h * nu**3 / c**2) * (h * nu / (k_B * T**2)) * (np.exp(h * nu / (k_B * T)) / (np.exp(h * nu / (k_B * T)) - 1)**2)

def rosseland_mean(x, nu, T):
    numerator = np.trapz(1 / x * dB_dT(nu, T), nu)
    denominator = np.trapz(dB_dT(nu, T), nu)
    return 1 / (numerator / denominator)

def plasma_freq(T):
    nu_P = 5.9281e-19 / h * (T/1e4)**(3/8) * np.exp(- 8.0109e-19 / (4 * k_B * T))
    return nu_P

def xi(nu, T):
    nu_P = plasma_freq(T)
    nu_int = nu[nu < nu_P]
    return np.trapz(B(nu_int, T), nu_int) / np.trapz(B(nu, T), nu)

data = np.loadtxt(args.input_file, skiprows=2, unpack=True)

energy = data[0]
sigma = data[1] * 1e6
nu = energy / h_ev
print('Plasma frequency:', plasma_freq(T) * h_ev, 'eV')
print('effective ratio:', 1 - xi(nu, T))
omega_in_hz = 2 * np.pi * nu
eps_imag = sigma / (omega_in_hz * epsilon_0)
eps_real = kkr_single(energy[1] - energy[0], eps_imag, cshift=args.cshift)

eps_mod = np.sqrt(eps_real**2 + eps_imag**2)
k = np.sqrt((eps_mod - eps_real) / 2)
n = np.sqrt((eps_mod + eps_real) / 2)

lambda_in_m = c / nu
wave_number_in_cm_inv = energy / (c * h_ev * 1e2)
alpha = 4 * np.pi * k / lambda_in_m

alpha_r = rosseland_mean(alpha, nu, T)
n_r = rosseland_mean(n, nu, T)
print('alpha_r:', alpha_r * 1e-2, 'cm^-1')
print('n_r:', n_r)
stefan_boltzmann = 5.67e-8
k_rad = 16 * n_r**2 * stefan_boltzmann * T**3 / (3 * alpha_r)
print('k_rad:', k_rad)

k_rad_nv = 16 * n**2 * stefan_boltzmann * T**3 / (3 * alpha)

fig, ax = mpt.init_plot_single()
ax.plot(energy, eps_real, label = r'Real part $\varepsilon_1$')
ax.plot(energy, eps_imag, label = r'Imaginary part $\varepsilon_2$')
mpt.legend(ax)
ax.minorticks_on()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Dielectric function')

mpt.savepdf(directory + 'epsilon')
np.savetxt(directory + 'epsilon.txt', np.column_stack((energy, eps_real, eps_imag)), fmt='%.6f')
# plt.show()
plt.cla()

ax.plot(energy, n, label = 'Real part (refractive index $n$)')
ax.plot(energy, k, label = 'Imaginary part (extinction coefficient $k$)')
mpt.legend(ax)
ax.minorticks_on()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Optical constants')

mpt.savepdf(directory + 'optical_constants')
np.savetxt(directory + 'optical_constants.txt', np.column_stack((energy, n, k)), fmt='%.6f')
# plt.show()
plt.close()

fig, axes = mpt.init_plot_double(sharey=False, wspace=0.1)
axes[0].plot(wave_number_in_cm_inv[5:100], alpha[5:100] * 1e-2)
axes[1].plot(lambda_in_m[5:50] * 1e9, alpha[5:50] * 1e-2)
axes[0].set_xlabel('Wave number (cm$^{-1}$)')
axes[0].set_yscale('log')
axes[1].set_xlabel('Wavelength (nm)')
axes[0].set_ylabel('Absorption coefficient (cm$^{-1}$)')
secax = axes[0].secondary_xaxis('top', functions=(lambda x: x * 299792458 * 4.135667696e-15 * 1e2, lambda x: x / (299792458 * 4.135667696e-15 * 1e2)))
secax.set_xlabel('Energy (eV)')
secax.minorticks_on()
for ax in axes:
    ax.minorticks_on()
    # ax.set_yscale('log')

mpt.savepdf(directory + 'alpha')
np.savetxt(directory + 'alpha.txt', np.column_stack((energy, alpha)), fmt='%.6f')
plt.close()

fig, ax = mpt.init_plot_single()
ax.plot(energy, k_rad_nv, label = 'Non-averaged')
ax.axhline(y=k_rad, color='r', linestyle='--', label = 'Averaged')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Radiative thermal conductivity (W/mK)')
mpt.legend(ax)
mpt.savepdf(directory + 'k_rad')
np.savetxt(directory + 'k_rad.txt', np.column_stack((energy, k_rad_nv)), fmt='%.6f')
np.savetxt(directory + 'k_rad_averaged.txt', np.array([k_rad]), fmt='%.6f')
# plt.show()
