#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

dirnames = ['nvt', 'nvt2', 'nvt3']
endsteps = np.loadtxt('ENDSTEPS', unpack=True, dtype=str)

dataout = open('full_nvt.dump', 'w')
for i in range(len(dirnames)):
    datain = open(dirnames[i]+'/nvt.dump', 'r')
    line = datain.readline()
    dataout.write(line)
    while line:
        line = datain.readline()
        if line.rstrip() == "ITEM: TIMESTEP":
            nextline = datain.readline()
            if nextline.rstrip() == endsteps[i]:
                break
            else:
                dataout.write(line)
                dataout.write(nextline)
        else:
            dataout.write(line)
    datain.close()
dataout.close()