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
 ./r_composition.py -r data/data.Rdata \
	-v EnvType -r1 Kingdom -s1 Bacteria -r2 Phylum -n 9 \
	-l test/phylo_compo.log -o test/phylo_compo.nb.html
