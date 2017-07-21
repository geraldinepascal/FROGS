#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
  rm -r test/*
else
  mkdir test
fi 


# python r_manova.py --rdata data/data.Rdata \
#                    --varExp EnvType \
#                    --distance-matrix data/Unifrac.tsv \
#                    --html test/phylo_manova.html \
#                    -l test/phylo_manova.log

python r_manova.py --rdata data/data.Rdata \
                   --varExp "EnvType + FoodType" \
                   --distance-matrix data/Unifrac.tsv \
                   --html test/phylo_manova.html \
                   -l test/phylo_manova.log
