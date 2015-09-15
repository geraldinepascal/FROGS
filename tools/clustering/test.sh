#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance.biom --output-fasta test/seeds.fasta --output-compo test/swarms_composition.tsv \
                --log-file test/log.txt

./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance2.biom --output-fasta test/seeds2.fasta --output-compo test/swarms_composition2.tsv \
                --log-file test/log2.txt \
                --denoising