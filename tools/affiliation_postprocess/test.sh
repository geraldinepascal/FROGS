#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
   rm test/*
fi

./affiliation_postprocess.py -b data/affiliations.biom -f data/filters.fasta -r data/Unite_extract_ITS1.fasta \
		--identity 65 --coverage 65 \
		--output-biom test/affiliations_postprocessed.biom --output-compo test/affiliations_postprocessed.compo \
		--output-fasta test/affiliations_postprocessed.fasta --log-file test/affiliations_postprocess.log --debug
