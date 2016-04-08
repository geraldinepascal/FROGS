#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# BIOM without affiliation
./biom_to_stdBiom.py \
 --input-biom data/abundance.biom \
 --output-biom test/abundance_1.biom \
 --output-metadata test/metadata_1.tsv \
 --log-file test/log_1.txt

# BIOM with affiliation
./biom_to_stdBiom.py \
 --input-biom data/affiliation.biom \
 --output-biom test/abundance_2.biom \
 --output-metadata test/metadata_2.tsv \
 --log-file test/log_2.txt
