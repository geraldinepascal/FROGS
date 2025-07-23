#!/bin/sh
# conda activate fogs@5.0.2
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

if [ -d test ]
then
	rm -r test/*
else
	mkdir test
fi 
 ./phyloseq_composition.py --phyloseq-rdata data/data.Rdata \
	--var-exp EnvType --taxa-rank-1 Kingdom --taxa-set-1 Bacteria --taxa-rank-2 Phylum --number-of-taxa 9 \
	--log-file test/phylo_compo.log --html test/phylo_compo.nb.html
