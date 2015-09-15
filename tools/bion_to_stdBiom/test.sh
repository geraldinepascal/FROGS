#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

./biom_to_stdBiom.py \
 --input-biom data/abundance.biom \
 --output-biom test/abundance_1.biom \
 --output-metadata test/metadata_1.tsv \
 --log-file test/log_1.txt

./biom_to_stdBiom.py \
 --input-biom data/affiliation.biom \
 --output-biom test/abundance_2.biom \
 --output-metadata test/metadata_2.tsv \
 --log-file test/log_2.txt