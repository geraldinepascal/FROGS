#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 

python r_clustering.py \
    --rdata data/data.Rdata \
    --varExp EnvType \
    --distance-matrix data/Unifrac.tsv \
    --html test/phylo_clustering.html \
    -l test/phylo_clustering.log
