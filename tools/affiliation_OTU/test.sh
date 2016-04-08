#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./affiliation_OTU.py --reference data/db.fasta \
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/affiliation.biom \
                     --summary test/summary.html --log-file test/aff.log
