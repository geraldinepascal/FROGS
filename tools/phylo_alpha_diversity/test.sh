#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

if [ -d test ]
then
	rm -r test/*
else
	mkdir test
fi 

python r_alpha_diversity.py -d data/data.Rdata \
	-v EnvType -m Observed Chao1 Shannon  \
	-l test/phylo_alpha.log -o test/phylo_alpha.html -a test/phylo_alpha.tsv
