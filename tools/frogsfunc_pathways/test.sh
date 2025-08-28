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
    --nb-cpus 2 \
    --input-tsv data/26-frogsfunc_functions_unstrat_abundance_EC.tsv \
    --normalisation \
    --strat-contrib \
    --input-asv-copy-norm data/26-frogsfunc_functions_marker_norm.tsv \
    --input-fun-copy data/EC_copynumbers_predicted.tsv \
    --output-pathways-abund test/27-frogsfunc_pathways_unstrat.tsv \
    --output-pathways-contrib test/27-frogsfunc_pathways_strat.tsv \
    --output-pathways-predictions test/27-frogsfunc_pathways_predictions.tsv \
    --output-pathways-abund-per-seq test/27-frogsfunc_pathways_unstrat_per_seq.tsv \
    --log-file test/frogsfunc_pathways_EC.log \
    --html test/27-frogsfunc_pathways_summary.html