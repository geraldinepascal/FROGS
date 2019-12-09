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

./r_alpha_diversity.py -r data/data.Rdata \
	-v EnvType -m Observed Chao1 Shannon  \
	-l test/phylo_alpha.log -o test/phylo_alpha.nb.html -a test/phylo_alpha.tsv --debug
