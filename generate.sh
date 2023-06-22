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

atomsk $1 -select random $2 Mg -substitute Mg He -properties types.txt tmp.lmp
atomsk tmp.lmp tmp.pdb

distance.py -i tmp.pdb.pdb
rm types.txt tmp.pdb.pdb

declare -a index_array=()

while IFS= read -r line
do
  index_array+=("$line")
done < "INDEX"

result=""

for (( i=0; i<${#index_array[@]}; i++ ))
do
    result+="${index_array[$i]}"
    result+=" 0 0 0.8 -add-atom H relative "
done

for (( i=0; i<${#index_array[@]}; i++ ))
do
    result+="${index_array[$i]}"
    if (( i != ${#index_array[@]}-1 ))
    then
        result+=" 0 0 -0.8 -add-atom H relative "
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

atomsk tmp.lmp -add-atom H relative $result 0 0 -0.8 -properties types-H.txt test.lmp
atomsk test.lmp -rmatom He brg-H.lmp
rm test.lmp
rm types-H.txt