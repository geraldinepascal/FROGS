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

./FPStep3.py \
 -b data/FPStep1.biom \
 -f data/FPStep2_predicted_functions.tsv \
 -m data/FPStep2_marker_nsti.tsv \
 --function-abund test/test_FPStep3_pred_abund_unstrat.tsv \
 --seqtab test/test_FPStep3_seqtab.tsv \
 --weighted test/test_FPStep3_weighted.tsv \
 -e test/test_FPStep3_excluded.txt \
 -l test/test_FPStep3.log \
 --html test/test_FPStep3.html