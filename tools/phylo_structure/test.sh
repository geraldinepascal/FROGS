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

./r_structure.py \
     --rdata data/data.Rdata \
     --varExp EnvType \
     --ordination-method MDS \
     --distance-matrix data/Unifrac.tsv \
     --html test/phylo_structure.nb.html \
     --log-file test/phylo_structure.log
