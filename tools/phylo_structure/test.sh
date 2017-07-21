#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 

python r_structure.py \
     --rdata data/data.Rdata \
     --varExp EnvType \
     --ordination-method MDS \
     --distance-matrix data/Unifrac.tsv \
     --html test/phylo_structure.html \
     --log-file test/phylo_structure.log
