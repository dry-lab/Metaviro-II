#!/bin/bash -e

source config/default_sh_values.cfg

while read -r model_file pattern port dest; do
    pkill -9 -f 'vw.*--port '$port
done < $model_pattern_port_file




