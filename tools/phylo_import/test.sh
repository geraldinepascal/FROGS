#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH


# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

# without normalisation and default ranks
mkdir -p test/
python r_import_data.py  \
	-n \
	-b data/chaillou.biom \
	-s data/sample_data.tsv \
	-t data/tree.nwk \
	-d test/phylo_import.Rdata -o test/phylo_import.html -l test/phylo_import.log
