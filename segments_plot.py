#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.integrate
from matplotlib import rcParams
import pickle as pkl
import argparse
import os

rcParams['font.size'] = '14'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

parser = argparse.ArgumentParser(description="Plot the time convergence of kappa")
parser.add_argument("--input_file", "-i", type=str, default="jh.dat",  help="log_lmp generated file")
parser.add_argument("--timestep", "-ts", type=float,  help=" timestep in fs, default 1fs")
parser.add_argument("--temperature", "-t", type=float,help='temperature in K')
parser.add_argument("--volume", "-v", type=float,help='volume in A3')
parser.add_argument("--configeration", "-conf", type=str, default='nve.lmp', help="read volume from this file")
parser.add_argument("--length", "-l", type=float, default=1, help='length of each segment')

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

kB   = 1.3806504e-23   # [J/K] Boltzmann
ev2j = 1.60218e-19
A2m  = 1.0e-10
ps2s = 1.0e-12
convert = ev2j*ev2j/ps2s/A2m
sample_rate = 1
scale = convert/kB/T/T/V*sample_rate*timestep
metal2SIps = convert/kB/T/T/V
file = args.input_file
J=np.loadtxt(file)
J = J*V
Jx, Jy, Jz = J[:,1], J[:,2], J[:,3]

def cumsum(xrange,Jx = Jx, Jy=Jy, Jz=Jz):
    JxJx = autocorr(Jx[xrange])
    print('Jx')
    JyJy = autocorr(Jy[xrange])
    print('Jy')
    JzJz = autocorr(Jz[xrange])
    print('Jz')
    JJ = (JxJx + JyJy + JzJz)/3
    cumsum_JxJx = scipy.integrate.cumtrapz(JxJx,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JyJy = scipy.integrate.cumtrapz(JyJy,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JzJz = scipy.integrate.cumtrapz(JzJz,initial=0)*scale; #np.insert(cumsum_JxJx, 0, 0)
    cumsum_JJ = (cumsum_JxJx + cumsum_JyJy + cumsum_JzJz)/3
    return cumsum_JJ, JJ*metal2SIps

l = args.length
num = len(Jx) // int(l * 1000 / timestep)

tcorr = np.array(range(int(l * 1000 / timestep)))*timestep # in ps
tcorr_idx = len(tcorr)
k_list, j_list, t_list = [], [], []
for i in range(num):
    print('Calculating round %d ...' % (i))
    if i == num - 1:
        k, j =cumsum(range(int(l * i * 1000 / timestep), len(Jx)))
    else:
        k, j =cumsum(range(int(l * i * 1000 / timestep), int(l * (1 + i) * 1000 / timestep)))
    k_list.append(k)
    j_list.append(j)
    t_list.append(len(k) * timestep / 1000)

with open('segments.pkl', 'wb') as file:
    data = [tcorr, j_list, k_list]
    pkl.dump(data, file)

fig,ax = plt.subplots(2,1,figsize=(8,12),sharex=True)
fig.subplots_adjust(hspace=0.05)

for index, j in enumerate(j_list):
    if index == len(j_list) - 1:
        rang = min(len(tcorr), len(j))
        ax[0].plot(tcorr[:rang], j[:rang], label="%d, %.2f ns" % (index + 1, t_list[index]))
    else:
        ax[0].plot(tcorr, j,label="%d, %.2f ns" % (index + 1, t_list[index]))
ax[0].set_ylabel('C(t)' + ' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1} \mathrm{ps}^{-1})$')
ax[0].plot(tcorr, np.ones(tcorr.shape)*0,'k--')
ax[0].set_xscale('log')
ax[0].set_ylim(bottom = - 0.1 * max(j), auto = True)

for index, k in enumerate(k_list):
    if index == len(k_list) - 1:
        rang = min(len(tcorr), len(k))
        ax[1].plot(tcorr[:rang], k[:rang], label="%d, %.2f ns" % (index + 1, t_list[index]))
    else:
        ax[1].plot(tcorr, k, label="%d, %.2f ns" % (index + 1, t_list[index]))
ax[1].set_xlabel('t (ps)')
ax[1].set_ylabel(r'$k $'+' '+ r'$(\mathrm{W} \mathrm{m}^{-1} \mathrm{K}^{-1})$')
ax[1].legend()
ax[1].set_xlim(min(tcorr[np.nonzero(tcorr)]), 0.5 * 1000 / 2)
ax[1].set_ylim(bottom = 0, top = 8, auto = True)
ax[1].set_xscale('log')
plt.savefig('segments.jpg', dpi=300, bbox_inches='tight')