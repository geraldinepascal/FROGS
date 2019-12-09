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

./r_manova.py --rdata data/data.Rdata \
                   --varExp "EnvType + FoodType" \
                   --distance-matrix data/Unifrac.tsv \
                   --html test/phylo_manova.nb.html \
                   -l test/phylo_manova.log
