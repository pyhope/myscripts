#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 11:02:36 2021
unit 
https://lammps.sandia.gov/doc/compute_heat_flux.html
directly analyze the log.properties file with Step, Jx, Jy, Jz 
@author: jiedeng
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 22:26:44 2021
analyze output of ave/correlate
@author: jiedeng
"""

#from util import reverse_readline
#import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.integrate
import pickle as pkl
from matplotlib import rcParams
rcParams['font.size'] = '14'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

parser = argparse.ArgumentParser(description="Plot contents from lammps log files")
parser.add_argument("--input_file", "-i",type=str, default="log.properties",  help="log_lmp generated file")
parser.add_argument("--timestep", "-ts",type=float,  help=" timestep in fs, default 1fs")
parser.add_argument("--temperature", "-t",type=float,help='temperature in K')
parser.add_argument("--volume", "-v",type=float,help='volume in A3')
parser.add_argument("-conf", "--configeration", type=str, default='conf.lmp', help="read volume from this file")
parser.add_argument("-s", "--store", default=False, action='store_true', help="Defualt:  Do not save data as outfile")
parser.add_argument("-of", "--outfile",type=str,default='log.kappa', help="out file name")
parser.add_argument("-sfig", "--store_fig", default=False, action='store_true', help="Defualt:  Do not save the figure as outfig")
parser.add_argument("-ofig", "--outfig",type=str,default='kappa.jpg', help="out figure name")

args = parser.parse_args()

def autocorr(a):
    b=np.concatenate((a,np.zeros(len(a))),axis=0)
    c= np.fft.ifft(np.fft.fft(b)*np.conjugate(np.fft.fft(b))).real
    d=c[:len(c)//2]
    d=d/(np.array(range(len(a)))+1)[::-1]
    return d

if args.temperature:
    T = args.temperature
else:
    infile = glob.glob('in*')[0]
    print('  Find ', infile)
    print("  ?? temperature not provided, parse from in file") 
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

sample_rate = 1       # for all mgsio3 syste, sample every single step

scale = convert/kB/T/T/V*sample_rate*timestep

metal2SIps = convert/kB/T/T/V

file = args.input_file

J=np.loadtxt(file) # heat current J is saved, but heat flux (energy x velocity) autocorrelation is saved in J0Jt, heat flux/V = heat flux
J = J*V
Jx, Jy, Jz = J[:,1], J[:,2], J[:,3]

#correlation time

def cumsum(xrange,Jx = Jx, Jy=Jy, Jz=Jz):
    JxJx = autocorr(Jx[xrange])
    JyJy = autocorr(Jy[xrange])
    JzJz = autocorr(Jz[xrange])
    JJ = (JxJx + JyJy + JzJz)/3
    cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3
    return cumsum_JJ, JJ*metal2SIps

tcorr = np.array(range(1000*1000))*timestep # in ps
tcorr_idx = range(len(tcorr))
    
fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)
fig.subplots_adjust(hspace=0.05)

k_t500, j_t500 =cumsum(range(1000*1000))
k_t1ns, j_t1ns =cumsum(range(2000*1000))
k_t2ns, j_t2ns =cumsum(range(4000*1000))
k_t4ns, j_t4ns =cumsum(range(len(Jx)))

ax[0].plot(tcorr, j_t500[tcorr_idx],label='0.5 ns')
ax[0].plot(tcorr, j_t1ns[tcorr_idx],label='1 ns')
ax[0].plot(tcorr, j_t2ns[tcorr_idx],label='2 ns')
ax[0].plot(tcorr, j_t4ns[tcorr_idx],label='4 ns')
#ax[0].set_ylabel('autocorrelation*V/kb/T^2  (W m-1 K-1 ps-1)')
ax[0].set_ylabel('C(t)' + ' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1} \mathrm{ps}^{-1})$')

#ax[0].grid(True)
ax[0].plot(tcorr, np.ones(tcorr.shape)*0,'k--')
ax[0].set_xscale('log')

ax[1].plot(tcorr, k_t500[tcorr_idx],label='0.5 ns')
ax[1].plot(tcorr, k_t1ns[tcorr_idx],label='1 ns')
ax[1].plot(tcorr, k_t2ns[tcorr_idx],label='2 ns')
ax[1].plot(tcorr, k_t4ns[tcorr_idx],label='4 ns')
#ax[1].set_xlabel('dt (ps)')
#ax[1].set_ylabel('thermal conductivity (W m-1 K-1)')
#ax[1].plot(tcorr, np.ones(tcorr.shape)*2.4,'k--')

ax[1].set_xlabel('t (ps)')
ax[1].set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')
#ax[1].grid(True)

ax[1].legend()
#ax[1].set_ylim([0,3.5])
#ax[1].set_xlim(min(tcorr),10)

ax[1].set_xscale('log')

if args.store_fig:
    plt.savefig(args.outfig, dpi=300, bbox_inches='tight')

if args.store:
    with open(args.outfile + '.pkl', 'wb') as file:
        data = [tcorr, j_t500, j_t1ns, j_t2ns, j_t4ns, k_t500, k_t1ns, k_t2ns, k_t4ns]
        pkl.dump(data, file)