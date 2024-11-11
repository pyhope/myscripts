#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 10:12:43 2022

submit new run with different volume, NSW, based on current old wrong

@author: jiedeng
"""
import os
from shutil import copyfile
import re
import fileinput
# import numpy as np
import argparse
from subprocess import call

parser = argparse.ArgumentParser()
parser.add_argument('--volume',"-v", type=float,help="new volume")
parser.add_argument('--target_path',"-tp", type=str,help="target path")

parser.add_argument('--nsw',"-n", type=int,help="new NSW, default=1")
parser.add_argument('--runvasp',"-rv", default=True, action='store_false', help='default, run vasp')
parser.add_argument('--subscript',"-s",default="sbatch ~/script/sub/sub.sh",help="command to submit vasp job, set qsubvasp as alias")
args   = parser.parse_args()

print("Build {0}".format(args.target_path))
os.mkdir(args.target_path)
copyfile('INCAR', os.path.join(args.target_path,'INCAR'))
copyfile('POTCAR', os.path.join(args.target_path,'POTCAR'))
copyfile('CONTCAR', os.path.join(args.target_path,'POSCAR'))
copyfile('KPOINTS', os.path.join(args.target_path,'KPOINTS'))
cwd = os.path.abspath(os.curdir)
os.chdir(args.target_path)

def key_word_replace(fname='INCAR',key_word_search='ENCUT',replacement_text='600'):
    """
    replace all the number for key word
    key_word_search : key word to search
    the problem come from 400.00 and 400 is different
    """
    tx_org = 'no original word found'
    tx_rep = 'no replacement done'    
    ## wild card is used here
#    for f_ele in os.listdir('.'):
#        if fnmatch.fnmatch(f_ele, fname):
#            f_full = f_ele
    f_full = os.path.abspath(fname)
    if os.path.isfile(f_full):
        with fileinput.FileInput(f_full, inplace=True, backup='-bak') as file:
            for line in file:
                if key_word_search in line:
                    
                    tx_org = line
                    #line=line.split()
                    #line_split = =line.split()
                    num = re.findall('\d*\.?\d+',line)
                    text_to_replace = ' '.join(map(str, num))
                    print(line.replace(text_to_replace, replacement_text))
                    tx_rep = replacement_text
                else:
                    print(line, end='')
#            os.rename(f_full+'-bak','old'+f_full)
    else:
        print('No'+ f_full+ 'found')
    return tx_org, tx_rep

if args.volume:
    print('modify the volume to {0}'.format(args.volume))
    call("python ~/script/format_converter.py -i POSCAR -s {0} -of vasp -o POSCAR".format(args.volume), shell=True)
if args.nsw:
    print('modify the NSW to {0}'.format(args.nsw))
    tx_org, tx_rep = key_word_replace(fname='INCAR',key_word_search='NSW', replacement_text=str(args.nsw))

if args.runvasp:
    call(args.subscript,shell=True)
    os.chdir(cwd)





