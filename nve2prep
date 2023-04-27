#!/bin/bash

temp=${1:-'4'}
press=${2:-'140'}

if test -e ./nvt
then
    echo ./nvt found!
else
    echo ./nvt not found!
    exit 1
fi

cd ./nvt

if test -e ./log.lammps
then
    echo log.lammps found!
else
    echo log.lammps not found!
    exit 1
fi

find_best_frame -st 10000 -t $temp -p $press

if test -e ./BEST_STEP
then
    echo Best step found!
else
    echo Best step not found!
    exit 1
fi

beststep=$(cat BEST_STEP)

mkdir ../nve2

extract_dump_frame -i nvt.dump -t $beststep
mv selected.dump ../nve2/

cd ..
sysname=`basename $PWD`
cp ../_inputs/nve/in.lammps nve2/
cp ../_inputs/nve/run.slurm nve2/

cd nve2
sed -i "s/job-name=nve/job-name=t${temp}p${press}-${sysname}/g" run.slurm
atomsk selected.dump nve.lmp
sed -i 's/1   1.00800000              # H/1   24.30500000             # Mg/g' nve.lmp
sed -i 's/2   4.00260200              # He/2   28.08500000             # Si/g' nve.lmp
sed -i 's/3   6.94000000              # Li/3   15.99900000             # O/g' nve.lmp
sed -i 's/4   9.01218200              # Be/4   1.00800000              # H/g' nve.lmp

#sbatch run.slurm
