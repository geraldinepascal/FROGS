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
    --strat-out \
    --marker-type 16S \
    --input-biom data/25-frogsfunc_placeseqs.biom \
    --input-fasta data/25-frogsfunc_placeseqs.fasta \
    --input-marker data/25-frogsfunc_marker.tsv \
    --input-tree data/25-frogsfunc_placeseqs_tree.nwk \
    --output-function-abund test/26-frogsfunc_functions_unstrat.tsv \
    --output-asv-norm test/26-frogsfunc_functions_marker_norm.tsv \
    --output-weighted test/26-frogsfunc_functions_weighted_nsti.tsv \
    --output-excluded test/26-frogsfunc_functions_excluded.txt \
    --output-contrib test/26-frogsfunc_functions_strat.tsv \
    --output-fasta test/26-frogsfunc_function.fasta \
    --output-biom test/26-frogsfunc_function.biom \
    --summary test/26-frogsfunc_functions_summary.html


