#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

./affiliation_OTU.py --reference data/SILVA_119_SSURef_Nr99_tax_silva.fasta \
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/affiliation.biom \
                     --summary test/summary.tsv --log-file test/aff.log
