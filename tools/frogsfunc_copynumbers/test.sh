#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./frogsfunc_copynumbers.py \
 -b data/25-frogsfunc_placeseqs.biom\
 -t data/25-frogsfunc_placeseqs_tree.nwk\
 -m test/26-frogsfunc_copynumbers_marker.tsv \
 -o test/26-frogsfunc_copynumbers_predicted_functions.tsv \
 --html test/26-frogsfunc_copynumbers_summary.html