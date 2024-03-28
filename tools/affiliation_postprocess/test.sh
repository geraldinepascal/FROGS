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

./affiliation_postprocess.py --input-biom data/affiliations.biom --input-fasta data/filters.fasta --reference data/Unite_extract_ITS1.fasta \
        --identity 65 --coverage 65 \
        --output-biom test/affiliations_postprocessed.biom --output-compo test/affiliations_postprocessed.compo \
        --output-fasta test/affiliations_postprocessed.fasta --log-file test/affiliations_postprocessed.log

# do not perform OTU aggregation based on taxonomies that include specific taxons
# mostly usefull with something like "unknown species"
./affiliation_postprocess.py --input-biom data/affiliations.biom --input-fasta data/filters.fasta --reference data/Unite_extract_ITS1.fasta \
        --identity 65 --coverage 65 --taxon-ignore "s__Cortinarius_papaver" \
        --output-biom test/affiliations_postprocessed_taxIgnore.biom --output-compo test/affiliations_postprocessed_taxIgnore.compo \
        --output-fasta test/affiliations_postprocessed_taxIgnore.fasta --log-file test/affiliations_postprocessed_taxIgnore.log
