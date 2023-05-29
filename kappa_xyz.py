#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from util import reverse_readline
import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.integrate
import pickle as pkl
from matplotlib import rcParams

rcParams['font.size'] = '16'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

parser = argparse.ArgumentParser(description="Plot the time convergence of kappa")
parser.add_argument("--input_file", "-i", type=str, default="jh.dat",  help="log_lmp generated file")
parser.add_argument("--timestep", "-ts", type=float, help=" timestep in fs, default 1fs")
parser.add_argument("--total_time", "-tt", type=float, default=100, help="total time for plotting in ps, default 100 ps")
parser.add_argument("--temperature", "-t", type=float,help='temperature in K')
parser.add_argument("--volume", "-v", type=float,help='volume in A3')
parser.add_argument("--configeration", "-conf", type=str, default='nve.lmp', help="read volume from this file")
parser.add_argument("--outfile", "-o", type=str,default='kappa', help="out file name")
parser.add_argument("--outfig", "-ofig", type=str,default='kappa.jpg', help="out figure name")

args = parser.parse_args()

def autocorr(a):
    b=np.concatenate((a,np.zeros(len(a))),axis=0)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)//2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d

if args.temperature:
    T = args.temperature
elif os.path.exists('Temp.txt'):
    print("  ?? temperature not provided, parse from Temp.txt") 
    with open("Temp.txt", "r") as file:
        T = float(file.read())
        print('  ** T = ', T)
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? temperature not provided, parse from in file. Warning: The results may be not accurate!!!!") 
    fp = open(infile)
    ins = fp.readlines()
    for line in ins:
        line=line.split('#')[0].split()
        if len(line)>0 and line[0] ==   'variable' and line[1] == 'T' and line[2] == 'equal':
            T = float(line[3])
            break
#            except:
#                print('No T found in ', infile)                   
    print('  ** T = ', T)

if args.timestep:
    timestep = args.timestep*1e-3#1e-3    # in ps this is different from post_corr
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? timestep not provided, parse from in file") 
    fp = open(infile)
    ins = fp.readlines()
    for line in ins:
        line=line.split('#')[0].split()
        if len(line)>0 and line[0] ==   'variable' and line[1] == 'dt' and line[2] == 'equal':
            timestep = float(line[3])
            break
#            except:
#                print('No T found in ', infile)                   
    print('  ** ts = ', timestep)
        
if args.volume:
    V = args.volume
else:
    infile = args.configeration
    fp = open(infile)
    ins = fp.readlines()
    print("  ?? vol not provided, parse from the configration file") 
    xhi = False
    yhi = False
    zhi = False
    for line in ins:
        if 'xlo' in line:
            xlo = float(line.split()[0])
            xhi = float(line.split()[1])
        if 'ylo' in line:
            ylo = float(line.split()[0])
            yhi = float(line.split()[1])
        if 'zlo' in line:
            zlo = float(line.split()[0])
            zhi = float(line.split()[1])
        if xhi and yhi and zhi:
            break
        
    V = (xhi - xlo)*(yhi - ylo)*(zhi - zlo)
    print(' ** volume = ', V)

# convert from LAMMPS real units to SI
kB   = 1.3806504e-23   # [J/K] Boltzmann
ev2j = 1.60218e-19
A2m  = 1.0e-10
ps2s = 1.0e-12
convert     = ev2j*ev2j/ps2s/A2m
scale = convert/kB/T/T/V*timestep
metal2SIps = convert/kB/T/T/V

file = args.input_file

J=np.loadtxt(file) # heat current J is saved, but heat flux (energy x velocity) autocorrelation is saved in J0Jt, heat flux/V = heat flux
J = J*V
Jx, Jy, Jz = J[:,1], J[:,2], J[:,3]
t = np.arange(len(Jx)) * timestep

#correlation time
tcorr = np.array(range(int(args.total_time/timestep)+1)) * timestep # in ps
tcorr_len = len(tcorr)

JxJx = autocorr(Jx)
print('Jx done!')
JyJy = autocorr(Jy)
print('Jy done!')
JzJz = autocorr(Jz)
print('Jz done!')

cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)

jxyz_full = np.array([JxJx, JyJy, JzJz]) * metal2SIps
kxyz_full = np.array([cumsum_JxJx, cumsum_JyJy, cumsum_JzJz])

jxyz = jxyz_full[:,:tcorr_len]
kxyz = kxyz_full[:,:tcorr_len]

j = np.mean(jxyz, axis=0)
k = np.mean(kxyz, axis=0)

fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)
fig.subplots_adjust(hspace=0.05)

ax[0].plot(tcorr, j, c='k', label='average (%d ns)' % (t[-1] / 1000))
ax[0].plot(tcorr, jxyz[0], c='C1', label='x')
ax[0].plot(tcorr, jxyz[1], c='C2', label='y')
ax[0].plot(tcorr, jxyz[2], c='C0', label='z')
ax[0].set_ylabel('C(t)' + ' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1} \mathrm{ps}^{-1})$')

ax[0].plot(tcorr, np.ones(tcorr.shape)*0,'k--')
ax[0].set_xscale('log')
ax[0].legend(fancybox=False, edgecolor='black')

ax[1].plot(tcorr, k, c='k', label='average (%.1f)' % (t[-1] / 1000))
ax[1].plot(tcorr, kxyz[0], c='C1', label='x')
ax[1].plot(tcorr, kxyz[1], c='C2', label='y')
ax[1].plot(tcorr, kxyz[2], c='C0', label='z')

ax[1].set_xlabel('t (ps)')
ax[1].set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')

ax[1].set_xscale('log')

plt.savefig(args.outfig, dpi=300, bbox_inches='tight')

with open(args.outfile + '.pkl', 'wb') as file:
    print('Total time for plotting:', tcorr[-1])
    data = [tcorr, jxyz, kxyz]
    pkl.dump(data, file)

with open(args.outfile + '_full.pkl', 'wb') as file:
    print('Total time for autocorrelation:', t[-1])
    data = [t, jxyz_full, kxyz_full]
    pkl.dump(data, file)