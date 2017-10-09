#!/bin/sh
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

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
	--rdata test/phylo_import.Rdata -o test/phylo_import.html -l test/phylo_import.log
