#!/bin/bash -e

usage="$(basename $0) [-h] [-f model_pattern_port_file] [-l model_location] [-n nb_categories] [-b hash_table_size]\n"
created_models=""

source config/default_sh_values.cfg

while getopts "hl:f:n:" opt; do
    case "$opt" in
	h) printf %b "$usage"; exit 0 ;;
	l) model_location=$OPTARG ;;
	f) model_pattern_port_file=$OPTARG ;;
	n) nb_cat=$OPTARG ;;
	b) vw_b_opt=$OPTARG ;;
    esac
done
shift $((OPTIND-1))

[ -z "${nb_cat##*[!0-9]*}" ] && error_exit '-n argument shoud be an integer'

test_existing_file $model_pattern_port_file
test_regular_file $model_pattern_port_file

while read -r model_file pattern port dest; do
    model=$model_location/$model_file

    if [ ! -e "$model" ]; then
	printf "*** File $model does not exist: generation of an empty model ***\n"
	init_model=$model.init

	touch $init_model
	for i in $(seq $nb_cat); do
	    printf "$i |\n" >> $init_model
	done && created_models="$created_models $init_model"

	vw -b $vw_b_opt -d $init_model -f $model -c -k --oaa $nb_cat --passes 20 --holdout_off
    fi

    test_regular_file $model && vw -i $model --daemon --num_children 1 --quiet --port $port
done < $model_pattern_port_file
