#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

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

# arbre pynast
mkdir -p test/pynast
python tree.py -i data/sequences.fasta \
	-t data/otus_pynast.fasta \
    -b data/sequences.biom \
    -s test/pynast/summary.html \
	-o test/pynast/frogs.nwk -l test/pynast/tree.log

# arbre mafft
mkdir -p test/mafft
python tree.py -i data/sequences.fasta \
    -b data/sequences.biom \
    -s test/mafft/summary.html \
	-o test/mafft/frogs.nwk -l test/mafft/tree.log

