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
    --nb-cpus 2 \
    --strat-contrib \
    --marker-type 16S \
    --functions EC COG \
    --input-biom data/25-frogsfunc_placeseqs.biom \
    --input-fasta data/25-frogsfunc_placeseqs.fasta \
    --input-marker-copy data/25-frogsfunc_marker_copy_per_asv.tsv \
    --input-tree data/25-frogsfunc_placeseqs_tree.nwk \
    --prefix-function-abund test/26-frogsfunc_functions_unstrat_abundance \
    --prefix-contrib  test/26-frogsfunc_functions_strat_contrib_and_abundance \
    --output-asv-copy-norm test/26-frogsfunc_functions_marker_norm.tsv \
    --output-weighted-nsti test/26-frogsfunc_functions_weighted_nsti.tsv \
    --output-excluded test/26-frogsfunc_functions_excluded.txt \
    --output-fasta test/26-frogsfunc_functions.fasta \
    --output-biom test/26-frogsfunc_functions.biom \
    --html test/26-frogsfunc_functions_summary.html \
    --log-file test/test_EC_COG.log --debug