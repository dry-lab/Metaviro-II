#!/bin/bash

path='vw/'
init_mod=${path}init_model.vw
nbcat=4
port=${1:-1234}

# create an empty model
echo -n > $init_mod
for i in $(seq $nbcat); do echo $i' |' >> $init_mod; done
vw -d ${path}train.vw -f ${path}init_model.vw -c -k --oaa $nbcat --passes 20 --holdout_off
# run vw daemon
vw -i ${path}init_model.vw --daemon --quiet --port $port -f ${path}model.vw
