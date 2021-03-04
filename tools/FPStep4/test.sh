#!/bin/sh

FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

# Dans mon ordi
# cd/Users/moussa/FROGS_moussa/tools/FPStep4/DonneesPourTest
python3 ../FPStep4.py -i pred_met1.tsv  --path_abund patth.tsv
