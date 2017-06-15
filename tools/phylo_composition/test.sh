#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
	rm -r test/*
else
	mkdir test
fi 

python2.7 r_composition.py -d data/data.Rdata \
	-v EnvType -r1 Kingdom -s1 Bacteria -r2 Phylum -n 5 \
	-l test/phylo_compo.log -o test/phylo_compo.html
