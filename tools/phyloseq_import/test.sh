#!/bin/sh
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$FROGS_DIR/app:$PATH
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
./phyloseq_import_data.py  \
	--normalisation \
	--biomfile data/chaillou.biom \
	--samplefile data/sample_data.tsv \
	--treefile data/tree.nwk \
	--rdata test/phylo_import.Rdata --html test/phylo_import.nb.html --log-file test/phylo_import.log

# on unstandardized biom
./phyloseq_import_data.py  \
	--normalisation \
	--biomfile data/frogs.biom \
	--samplefile data/frogs_sample_data.tsv \
	--rdata test/frogs_import.Rdata --html test/frogs_import.nb.html --log-file test/frogs_import.log
