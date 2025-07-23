#!/bin/sh
# conda activate frogs@5.0.2
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 

./phyloseq_clustering.py \
    --phyloseq-rdata data/data.Rdata \
    --var-exp EnvType \
    --beta-distance-matrix data/Unifrac.tsv \
    --html test/phylo_clustering.nb.html \
    --log-file test/phylo_clustering.log
