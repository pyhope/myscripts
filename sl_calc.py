#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

fes=np.genfromtxt("fes.dat")
N=fes.shape[0]
freeEnergyLiquid=np.amin(fes[:,1][:int(N/2)])
freeEnergySolid=np.amin(fes[:,1][int(N/2):])
deltaFreeEnergy=freeEnergySolid-freeEnergyLiquid

integration = np.sum(fes[:,1][:int(N/2)]) - np.sum(fes[:,1][int(N/2):])
print("%.2f" % (deltaFreeEnergy) + " %.2f" % (integration), end=' ')

if deltaFreeEnergy < 0:
    print("solid")
else:
    print("liquid")