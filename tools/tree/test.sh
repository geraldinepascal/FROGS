#!/bin/sh
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test
fi

##############
# FROGS_tree #
##############

# # arbre mafft
mkdir -p test/mafft
./tree.py -i data/sequences.fasta \
    -b data/sequences.biom \
    -s test/mafft/summary.html \
    -o test/mafft/frogs.nwk -l test/mafft/tree.log
