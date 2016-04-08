#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# Input BIOM
./remove_chimera.py \
                    --input-fasta data/swarm_seed_sequences.fasta --input-biom data/swarm_abundance.biom \
                    --non-chimera test/non_chimera.fasta --out-abundance test/abundance.biom --summary test/summary.html --log-file test/log.txt \
                    --nb-cpus 2

# Input count
./remove_chimera.py \
                    --input-fasta data/swarm_seed_sequences.fasta --input-count data/seq_count.tsv \
                    --non-chimera test/non_chimera2.fasta --out-abundance test/abundance2.tsv --summary test/summary2.html --log-file test/log2.txt \
                    --nb-cpus 2
