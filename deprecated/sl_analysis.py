#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import my_pyplot as mpt

COLVAR=np.genfromtxt("COLVAR")

fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(COLVAR[:,0],COLVAR[:,2])
ax.set_xlabel("Time (ps)")
ax.set_ylabel("# of crystal-like atoms")
mpt.savejpg('fig1')

plt.cla()
fes=np.genfromtxt("fes.dat")
# fes=np.genfromtxt("fes/fes_20.dat")
ax.plot(fes[:,0],fes[:,1]-np.amin(fes[:,1]))
ax.set_xlabel("Collective variable")
ax.set_ylabel("Free energy (kJ/mol)")
mpt.savejpg('fig2')

N=fes.shape[0]
freeEnergyLiquid=np.amin(fes[:,1][:int(N/2)])
freeEnergySolid=np.amin(fes[:,1][int(N/2):])
deltaFreeEnergy=freeEnergySolid-freeEnergyLiquid
#print("Free energy difference: " + "%.2f" % (deltaFreeEnergy) + " kJ/mol")
#print("Free energy difference: " + "%.2f" % (deltaFreeEnergy) + " kJ/mol or " + "%.2f" % (deltaFreeEnergy/14.134) + " kT.")

if deltaFreeEnergy < 0:
    print("%.2f" % (deltaFreeEnergy) + " solid")
else:
    print("%.2f" % (deltaFreeEnergy) + " liquid")