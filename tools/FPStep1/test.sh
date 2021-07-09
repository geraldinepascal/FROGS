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

python2.7 placement.py -s data/test.fasta -o test/out5.tree -r data/ref_tree_picrust2/default_files/prokaryotic/pro_ref -b test.biom1 --min_align 0.8
# Dans mon ordi
#python3 ../../FPStep1.py -s test_16S.fasta -o out5.tree -r /Users/moussa/FROGS_moussa/tools/FPStep1/data/ref_tree_picrust2/default_files/prokaryotic/pro_ref -b data/chemerin_16S/table.biom --min_align 0.8

