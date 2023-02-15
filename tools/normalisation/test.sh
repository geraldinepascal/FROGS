#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

#mkdir -p test/test_N1000
#./normalisation.py -n 1000 -i data/test_3OTU_2sample_change.biom -b test/test_N1000/test_N1000.biom \
#                         -s test/test_N1000/summary.html -l test/test_N1000/log.txt 

#mkdir -p test/test_M
#./normalisation.py -m -i data/test_3OTU_2sample_change.biom -b test/test_M/test_M.biom \
#                         -s test/test_M/summary.html -l test/test_M/log.txt

#mkdir -p test/test_aff_M
#./normalisation.py -m -i data/seed_abondance_affiliation2.biom -b test/test_aff_M/test_aff_M.biom \
#                         -s test/test_aff_M/summary.html -l test/test_aff_M/log.txt

./normalisation.py \
 -n 1000 \
 -i data/non_chimera_500sp.biom \
 -f data/non_chimera_500sp.fasta \
 -b test/test_aff_N1000.biom \
 -o test/test_aff_N1000.fasta \
 -s test/summary.html \
 -l test/log.txt