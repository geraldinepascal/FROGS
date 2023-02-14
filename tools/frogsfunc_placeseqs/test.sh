#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./frogsfunc_placeseqs.py \
 -i data/frogsfunc.fasta \
 -b data/frogsfunc.biom \
 -p sepp \
 -o test/25-frogsfunc_placeseqs_tree.nwk \
 -e test/25-frogsfunc_placeseqs_excluded.txt \
 -s test/25-frogsfunc_placeseqs.fasta \
 -m test/25-frogsfunc_placeseqs.biom\
 -c test/25-frogsfunc_placeseqs_closests_ref_sequences.txt \
 -om test/25-frogsfunc_marker.tsv \
 --summary test/25-frogsfunc_placeseqs_summary.html \
 --debug
