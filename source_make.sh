#! /bin/bash

IRad=(  "0.9"  "1.16"  "1.52" "1.66")
Label=( "Li"  "Na"  "K"  "Rb" )
RO="1.52"

for ((n=0; n<4; n++)); do
    L=${Label[$n]}
    R=$(echo  "$RO + ${IRad[$n]}" | bc -l)
    rm -r -f $L
    mkdir $L
    for ((i=4; i<=14; i++)); do
        awk '{printf " %15.10f", $1; if(NR % 3 ==0){printf "\n";}}' pack.3.$i.txt | sed "s/^/$L  /g" | sed "s/\$/   $R/g" > $L/pack_3.$i.vecfle
    done
done


