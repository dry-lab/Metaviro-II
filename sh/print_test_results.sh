#!/bin/bash

# find ../log/models/${exp}/${pattern}_* -cnewer ../log/models/${exp}/${pattern}_21March_08\:00.mod \
# $(ls -tr ../log/models/$exp/${pattern}_* | awk -F "/" '{print $NF}' | awk -F "." '{print $1}')

cats=(archea bact euk virus)

for exp in empty ref
do
    echo '* '"$exp"
    for pattern in 111 110110 10000111
    do
	echo '** '"$pattern"
	for file in $(ls $1/${exp}/${pattern}_* | awk -F "/" '{print $NF}' | awk -F "." '{print $1}') 
	do
	    fdate=$(echo $file | awk -F "_" '{print substr($2,0,5)}')
	    fhour=$(echo $file | awk -F "_" '{print $3}' | awk -F '.' '{print $1}')
	    echo -ne '|'"${fdate}_$fhour"'|'
	    for ((i = 1; i < 5; i++))
	    do
		fscore=$(grep -ce "^$i" $2/$exp/${cats[i-1]}/${file}*)
		echo -ne "\t$fscore"'|'
	    done
	    echo
	done
    done
done


