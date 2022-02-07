#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH
export GENE_HIERARCHY_FILE=$FROGS_DIR/default_files/gene_family_hierarchy.tsv

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./frogsfunc_functions.py \
 -b data/25-frogsfunc_placeseqs.biom\
 -f data/26-frogsfunc_copynumbers_predicted_functions.tsv \
 -m data/26-frogsfunc_copynumbers_marker.tsv \
 --function-abund test/27-frogsfunc_functions_unstrat.tsv \
 --seqtab test/27-frogsfunc_functions_marker_norm.tsv \
 --weighted test/27-frogsfunc_functions_weighted_nsti.tsv \
 -e test/27-frogsfunc_functions_excluded.txt \
 --html test/27-frogsfunc_functions_summary.html