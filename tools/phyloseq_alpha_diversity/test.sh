#!/bin/sh
#conda activate frogs@5.0.2
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

if [ -d test ]
then
	rm -r test/*
else
	mkdir test
fi 

./phyloseq_alpha_diversity.py --phyloseq-rdata data/data.Rdata \
	--var-exp EnvType --alpha-measures Observed Chao1 Shannon  \
	--log-file test/phylo_alpha.log --html test/phylo_alpha.nb.html --output-alpha-tsv test/phylo_alpha.tsv --debug
