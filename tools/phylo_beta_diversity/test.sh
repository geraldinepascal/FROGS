#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 

python r_beta_diversity.py \
           --data data/data.Rdata \
           --varExp EnvType \
           --distance-methods bray,unifrac,wunifrac,euclidean \
           --output-dir test \
           --html test/phylo_beta.html \
           --log-file test/phylo_beta.log
