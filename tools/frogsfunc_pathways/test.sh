#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH
export DESCRIPTION_FILE=$FROGS_DIR/frogsfunc_suppdata/pathways_description_file.txt.gz
export PATHWAYS_HIERARCHY_FILE=$FROGS_DIR/frogsfunc_suppdata/pathways_hierarchy.tsv

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./frogsfunc_pathways.py \
 -i data/27-frogsfunc_functions_unstrat.tsv \
 -o test/28-frogsfunc_pathways_unstrat.tsv \
 -t test/28-frogsfunc_pathways_summary.html