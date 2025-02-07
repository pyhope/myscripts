#!/bin/bash

set -e

WDIR=`pwd`

Number_of_bands=$(awk 'NR==6 {print $3}' EIGENVAL)
Nuber_of_bands_minus_1=$(($Number_of_bands-1))

Fermi_energy=$(awk 'NR==6 {print $4}' DOSCAR)

echo "Number of bands: $Number_of_bands"
echo "Fermi energy: $Fermi_energy"
grep -A$Nuber_of_bands_minus_1 "^    1" EIGENVAL | grep -v "^--$" | awk -v Fermi_energy="$Fermi_energy" '{print $1, $2-Fermi_energy, $4*100.0, $5*100.0}'
