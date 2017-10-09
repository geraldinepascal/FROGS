#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# Clustering without denoising
./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance.biom --output-fasta test/seeds.fasta --output-compo test/swarms_composition.tsv \
                --log-file test/log.txt

# Clustering with denoising step
./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance2.biom --output-fasta test/seeds2.fasta --output-compo test/swarms_composition2.tsv \
                --log-file test/log2.txt \
                --denoising
