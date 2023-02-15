#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH
export GENE_HIERARCHY_FILE=$FROGS_DIR/frogsfunc_suppdata/gene_family_hierarchy.tsv

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./frogsfunc_functions.py \
 -b data/25-frogsfunc_placeseqs.biom \
 -i data/25-frogsfunc_placeseqs.fasta \
 -m data/25-frogsfunc_marker.tsv \
 -t data/25-frogsfunc_placeseqs_tree.nwk \
 --output-dir test \
 --marker-type 16S \
 --output-function 26-frogsfunc_copynumbers_functions.tsv \
 --output-function-abund 26-frogsfunc_functions_unstrat.tsv \
 --output-otu-norm 26-frogsfunc_functions_marker_norm.tsv \
 --output-biom 26-frogsfunc_function.biom \
 --output-fasta 26-frogsfunc_function.fasta \
 --output-weighted 26-frogsfunc_functions_weighted_nsti.tsv \
 --output-excluded 26-frogsfunc_functions_excluded.txt \
 --summary 26-frogsfunc_functions_summary.html --debug 

