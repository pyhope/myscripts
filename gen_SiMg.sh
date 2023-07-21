#!/bin/bash

set -e

if [ -z $2 ]; then
    echo Missing parameters!
    exit 1
fi

cat > types.txt << EOF
types
Mg  1
Si  2
O   3
He  4
EOF

atomsk $1 -select random $2 Si -substitute Si He -properties types.txt tmp.lmp
atomsk tmp.lmp tmp.pdb

distance.py -i tmp.pdb.pdb
rm types.txt tmp.pdb.pdb

declare -a index_array=()

while IFS= read -r line
do
  index_array+=("$line")
done < "INDEX"

result1="-add-atom H relative "

for (( i=0; i<${#index_array[@]}; i++ ))
do
    result1+="${index_array[$i]}"
    result1+=" 0.474 -0.101 0.654 -add-atom H relative "
done

for (( i=0; i<${#index_array[@]}; i++ ))
do
    result1+="${index_array[$i]}"
    if (( i != ${#index_array[@]}-1 ))
    then
        result1+=" -0.757 0.190 0.242 -add-atom H relative "
    else
        result1+=" -0.757 0.190 0.242"
    fi
done

result2="-add-atom Mg relative "
for (( i=0; i<${#index_array[@]}; i++ ))
do
    result2+="${index_array[$i]}"
    if (( i != ${#index_array[@]}-1 ))
    then
        result2+=" 0.101 -0.003 -0.320 -add-atom Mg relative "
    else
        result2+=" 0.101 -0.003 -0.320"
    fi
done

cat > types-H.txt << EOF
types
Mg  1
Si  2
O   3
H   4
He  5
EOF

atomsk tmp.lmp $result1 $result2 -properties types-H.txt test.lmp
atomsk test.lmp -rmatom He brg-H.lmp
rm test.lmp
rm types-H.txt