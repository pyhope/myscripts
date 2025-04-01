#!/bin/bash

unit=${1:-'A'}
input_filename=${2:-'md'}
output_filename=${3:-'traj'}

atomsk --unfold $input_filename.out exyz
ls *.xyz | sort -V > xyz.lst

if [ "$unit" == "B" ]; then
for i in $(cat xyz.lst); do
atomsk $i -unit Bohr Angstroms $i-tmp.exyz
mv $i-tmp.xyz $i
done
fi

atomsk --gather xyz.lst $output_filename.exyz
xargs rm < xyz.lst
rm xyz.lst
