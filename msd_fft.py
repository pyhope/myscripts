#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:25:43 2020
## msd direct is wrong with unwrapping
@author: jiedeng
"""
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str,help="input file")
parser.add_argument("--eles","-e",type = str,default = None, help="elements to analyze")
parser.add_argument("--format","-ft",type = str,default = 'LAMMPSDUMP', help="file format, e.g., LAMMPSDUMP, PDB")
parser.add_argument("--timestep","-ts",type = int,default = 1, help="timestep")

args   = parser.parse_args()


#def unwrap(idx):
#    #Bulow et al., 2020 unwrap method eqn3
#    newcoords = []
#    for i in range(len(frames)):
#        if i == 0:
#            newcoords.append(frames[i].positions[idx])
#        else:
#            curr = frames[i].positions[idx]
#            
#            newcoords.append(curr - np.floor((curr - newcoords[i-1])/box +.5)*box)
#
#    return np.array(newcoords)


def unwrap(idx):
    #Bulow et al., 2020 unwrap method, eqn1
    newcoords = []
    for i in range(len(frames)):
        if i == 0:
            newcoords.append(frames[i].positions[idx])
        else:
            wr_i = frames[i].positions[idx]
            wr_i1 = frames[i-1].positions[idx]
            un_i1 = newcoords[i-1]
            tmp =   np.floor((wr_i - wr_i1)/box +.5)*box
            newcoords.append(un_i1 + wr_i - wr_i1 - tmp)
            
    return np.array(newcoords)


def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_fft(r):
  N=len(r)
  D=np.square(r).sum(axis=1) 
  D=np.append(D,0) 
  S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
  Q=2*D.sum()
  S1=np.zeros(N)
  for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
  return S1-2*S2

def msd_straight_forward(r):
    shifts = np.arange(len(r))
    msds = np.zeros(shifts.size)    

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()

    return msds
if args.file:
    file = args.file
else:
    print("No file supplied, search for *dump file")
    files  = glob.glob("*dump")
    file = files[0]
    print("Find {0}; analyze {1}".format(files,file))

#file = '/Users/jiedeng/Documents/ml/deepmd-kit/my_example/6k/tmp/r1/20/mgsio3.dump'
u_md = mda.Universe(file,format=args.format)
u_md.transfer_to_memory()
frames = u_md.trajectory
box    = frames[0].dimensions[0:3]

eles = []
for atom in u_md.atoms.types:
    if not (atom in eles):
        eles.append(atom)

# find unique elements
msds = []
if not (args.eles):
    ele_sel = eles
else:
    ele_sel = args.eles.split('-')
    
for ele in ele_sel:
    idx = u_md.select_atoms('type {0}'.format(ele)).indices
    print(ele)
    tmp = np.zeros(len(u_md.trajectory))
    for i in idx:
        r    = unwrap(i)
        tmp += msd_fft(r)
#        tmp += msd_straight_forward(r)
    msds.append(tmp/len(idx))
print(msds)
print("start ")
np.savetxt("msd_fft.txt", np.array(msds).T, header='    '.join(ele_sel), fmt = '%2.6f')
print("end saving ")

plt.figure(dpi=300)
for i in range(len(ele_sel)):
    plt.loglog((np.array(list(range(len(msds[i]))))*args.timestep)[1:],msds[i][1:],label=eles[i])
#plt.legend()
#plt.ylim([1e-1,1e2])
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd_fft.png",bbox_inches='tight')

plt.close()
plt.figure(dpi=300)
for i in range(len(ele_sel)):
    plt.plot((np.array(list(range(len(msds[i]))))*args.timestep)[1:],msds[i][1:], c='k' ,label=eles[i])
#plt.legend()
#plt.ylim([1e-1,1e2])
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd_linear.png",bbox_inches='tight')
