#!/bin/bash

cats=(archea bact euk virus)

#for exp in empty ref
#do
for pattern in 1111111111 11011011011011 10001110010110111 1010110000111110001 11110100011100011 #111 110110 10000111
do
    for file in $(ls $1/${exp}/${pattern}_*) #find $1/${exp}/${pattern}_* -cnewer $1/${exp}/${pattern}_21March_08\:00.mod)
    do
	for ((i = 1; i < 5; i++))
	do
	    outputfile=$(echo $file | awk -F "/" '{print $NF}' | awk -F '.' '{print $1}').predict
	    grep -E "^$i" all_vw_$pattern \
		| vw -t -i $file \
		-p $2/${exp}/${cats[i-1]}/$outputfile
	done
    done
done
#done

