#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
start_time = time.time()
formatted_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
print('Program starts at: ', formatted_time)
print("Loading modules ...")
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import argparse
import glob
import multiprocessing
import os
module_time = time.time()
print("End after %.2f s" % (module_time - start_time))

parser = argparse.ArgumentParser()
parser.add_argument("--file","-f",type=str,help="input file")
parser.add_argument("--eles","-e",type = str, default = None, help="elements to analyze (split with -)")
parser.add_argument("--boundary","-b",type=str, default = None, help="boundaries to analyze (format: x-y-z in Angstrom)")
parser.add_argument("--format","-ft",type = str,default = 'LAMMPSDUMP', help="file format, e.g., LAMMPSDUMP, PDB")
parser.add_argument("--timestep","-ts",type = int,default = 1, help="timestep")
parser.add_argument("--nptfile","-nptf",type=str, default = './npt.lmp', help="npt.lmp file")
parser.add_argument("--npt_atom_style","-nptas",type=str, default = 'id type q x y z', help="npt.lmp file atom style")
parser.add_argument("--boundary2", "-b2", default=False, action='store_true', help="Defualt: calculate the center boundary")
parser.add_argument("--boundary3", "-b3", default=False, action='store_true', help="Defualt: calculate the center boundary")
parser.add_argument("--boundary4", "-b4", default=False, action='store_true', help="Defualt: calculate the center boundary")
parser.add_argument("--start","-s",type = int, default = 0, help="start frame")

args   = parser.parse_args()

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
    msdsx = np.zeros(shifts.size)
    msdsy = np.zeros(shifts.size)
    msdsz = np.zeros(shifts.size)

    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        diffx, diffy, diffz = diffs[:,0], diffs[:,1], diffs[:,2]
        sqdistx, sqdisty, sqdistz = np.square(diffx), np.square(diffy), np.square(diffz)
        msdsx[i], msdsy[i], msdsz[i] = sqdistx.mean(), sqdisty.mean(), sqdistz.mean()

    return msdsx, msdsy, msdsz

def task(i):
    r = unwrap(i)
    data_x = msd_straight_forward(r)[0]
    data_y = msd_straight_forward(r)[1]
    data_z = msd_straight_forward(r)[2]
    print(i, end=' ')
    return data_x, data_y, data_z

if args.file:
    file = args.file
else:
    print("No file supplied, search for *dump file")
    files  = glob.glob("*dump")
    file = files[0]
    print("Find {0}; analyze {1}".format(files,file))

print("Loading dump file ...")
u_md = mda.Universe(file, format=args.format)
loading_time = time.time()
print("End after %.2f s" % (loading_time - module_time))

print("Transferring to memory ...")
u_md.transfer_to_memory()
transfer_time = time.time()
print("End after %.2f s" % (transfer_time - loading_time))

frames = u_md.trajectory[args.start:]
box    = frames[0].dimensions[0:3]

print("Calculating MSD ...")
cores = int(os.getenv('SLURM_CPUS_PER_TASK'))
print('Number of CPU cores', cores)
msdx, msdy, msdz = [], [], []
eles = []
for atom in u_md.atoms.types:
    if not (atom in eles):
        eles.append(atom)
# find unique elements
if not (args.eles):
    ele_sel = eles
else:
    ele_sel = args.eles.split('-')
for ele in ele_sel:
    if not args.boundary:
        idx = u_md.select_atoms('type %s' % (ele)).indices
    else:
        bx, by, bz = args.boundary.split('-')
        u_npt = mda.Universe(args.nptfile, format='DATA', atom_style=args.npt_atom_style)
        box_dim = u_npt.dimensions[:3]
        bx2, by2, bz2 = str(box_dim[0] - float(bx)), str(box_dim[1] - float(by)), str(box_dim[2] - float(bz))
        if args.boundary2:
            by_new = str(box_dim[1] / 2 - float(by))
            by2_new = str(box_dim[1] - float(by_new))
            idx = u_npt.select_atoms('type '+ele+' and prop x > '+bx+' and prop x < '+bx2+' and (prop y < '+by_new+' or prop y > '+by2_new+') and prop z > '+bz+' and prop z < '+bz2).indices
        elif args.boundary3:
            by3 = str(box_dim[1]/2 - float(by))
            idx = u_npt.select_atoms('type '+ele+' and prop y > '+by+' and prop y < '+by3).indices
        elif args.boundary4:
            by4 = str(box_dim[1]/2 + float(by))
            idx = u_npt.select_atoms('type '+ele+' and prop y > '+by4+' and prop y < '+by2).indices
        else:
            idx = u_npt.select_atoms('type '+ele+' and prop x > '+bx+' and prop x < '+bx2+' and prop y > '+by+' and prop y < '+by2+' and prop z > '+bz+' and prop z < '+bz2).indices
    print()
    print('Elements:', ele)
    print('Number of atoms:',len(idx))
    tmpx = np.zeros(len(frames))
    tmpy = np.zeros(len(frames))
    tmpz = np.zeros(len(frames))
    with multiprocessing.Pool(cores) as pool:
        results = pool.map(task, idx)
    for result in results:
        tmpx += result[0]
        tmpy += result[1]
        tmpz += result[2]
    msdx.append(tmpx/len(idx))
    msdy.append(tmpy/len(idx))
    msdz.append(tmpz/len(idx))

    calc_time = time.time()
    print()
    print("End after %.2f s" % (calc_time - transfer_time))
    transfer_time = calc_time

print()
print("Saving to file ...")
x_arr, y_arr, z_arr = np.array(msdx).T, np.array(msdy).T, np.array(msdz).T
np.savetxt("msd_x.txt", x_arr, header='    '.join(ele_sel), fmt = '%2.6f')
np.savetxt("msd_y.txt", y_arr, header='    '.join(ele_sel), fmt = '%2.6f')
np.savetxt("msd_z.txt", z_arr, header='    '.join(ele_sel), fmt = '%2.6f')
np.savetxt("msd_fft.txt", x_arr + y_arr + z_arr, header='    '.join(ele_sel), fmt = '%2.6f')
sf_time = time.time()
print("End after %.2f s" % (sf_time - calc_time))

print("Plotting ...")
plt.figure(dpi=300)
for i in range(len(ele_sel)):
    x = (np.array(list(range(len(msdx[i]))))*args.timestep)[1:]
    plt.loglog(x,msdx[i][1:],label=eles[i]+', x')
    plt.loglog(x,msdy[i][1:],label=eles[i]+', y')
    plt.loglog(x,msdz[i][1:],label=eles[i]+', z')
plt.legend()
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd_xyz_log.png",bbox_inches='tight')
plt.cla()
for i in range(len(ele_sel)):
    msd_bulk = msdx[i] + msdy[i] + msdz[i]
    plt.loglog(x,msd_bulk[1:],label=eles[i], c='k')
plt.legend()
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd_log.png",bbox_inches='tight')

plt.close()
plt.figure(dpi=300)
for i in range(len(ele_sel)):
    plt.plot(x,msdx[i][1:], label=eles[i]+', x')
    plt.plot(x,msdy[i][1:], label=eles[i]+', y')
    plt.plot(x,msdz[i][1:], label=eles[i]+', z')
plt.legend()
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd_xyz.png",bbox_inches='tight')
plt.cla()
for i in range(len(ele_sel)):
    msd_bulk = msdx[i] + msdy[i] + msdz[i]
    plt.plot(x,msd_bulk[1:],label=eles[i], c='k')
plt.legend()
plt.xlabel('time (fs)')
plt.ylabel('MSD (A$^2$)')
plt.savefig("msd.png",bbox_inches='tight')

plot_time = time.time()
print("End after %.2f s" % (plot_time - sf_time))
print("Total time: %.2f s" % (plot_time - start_time))
