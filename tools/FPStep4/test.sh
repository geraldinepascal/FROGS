#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH
export GENE_HIERARCHY_FILE=$FROGS_DIR/default_files/gene_family_hierarchy.tsv
export DESCRIPTION_FILE=$FROGS_DIR/default_files/pathways_description_file.txt.gz
export PATHWAYS_HIERARCHY_FILE=$FROGS_DIR/default_files/pathways_hierarchy.tsv

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./FPStep4.py \
 -i data/FPstep3_pred_abund.tsv \
 -o test/test_FPStep4_pathways_abund.tsv \
 -l test/test_FPStep4.log \
 -t test/test_FPStep4.html