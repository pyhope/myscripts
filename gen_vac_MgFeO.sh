#!/bin/bash

set -e

if [ -z $3 ]; then
    echo Missing parameters!
    exit 1
fi

atomsk $1 -select above 111 y -select intersect below 121 y -select among random $2 Mg -rmatoms select tmp.lmp
atomsk tmp.lmp -select above 111 y -select intersect below 121 y -select among random $2 O -rmatoms select removed.lmp

cat > types.txt << EOF
types
Mg  1
Fe  2
O   3
EOF

atomsk removed.lmp -select random 10% Mg -substitute Mg Fe -properties types.txt MgFeO.lmp

rm tmp.lmp removed.lmp types.txt