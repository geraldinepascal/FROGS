#!/bin/sh
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 

./r_beta_diversity.py \
           --rdata data/data.Rdata \
           --varExp EnvType \
           --distance-methods bray,unifrac,euclidean \
           --matrix-outdir test \
           --html test/phylo_beta.nb.html \
           --log-file test/phylo_beta.log
