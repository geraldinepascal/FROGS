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

./phyloseq_alpha_diversity.py --rdata data/data.Rdata \
	--varExp EnvType --alpha-measures Observed Chao1 Shannon  \
	--log-file test/phylo_alpha.log --html test/phylo_alpha.nb.html --alpha-out test/phylo_alpha.tsv --debug
