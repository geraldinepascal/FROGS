#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./FPStep2.py \
 -b data/FPStep1.biom \
 -t data/FPStep1.tree \
 -m test/test_FPStep2_marker_nsti.tsv \
 -o test/test_FPStep2_predicted_functions.tsv \
 -l test/test_FPStep2.log \
 --html test/test_FPStep2.html