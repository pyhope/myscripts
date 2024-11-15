#!/bin/bash

filename=${1:-'traj'}

atomsk --unfold ppv.out exyz
ls *.xyz > xyz.lst
atomsk --gather xyz.lst $filename.exyz
xargs rm < xyz.lst
rm xyz.lst
