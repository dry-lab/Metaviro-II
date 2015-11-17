#!/bin/bash

path='../vw/'
nbcat=4
port=${1:-1234}

# create an empty model
vw -d ${path}train.vw -f ${path}init_model.vw -c -k --oaa $nbcat --passes 20 --holdout_off
# run vw daemon
vw -i ${path}init_model.vw --daemon --quiet --port $port -f ${path}model.vw
