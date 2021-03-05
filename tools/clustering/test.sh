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

# Clustering without denoising
echo Clustering without denoising
./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance.biom --output-fasta test/seeds.fasta --output-compo test/swarms_composition.tsv \
                --log-file test/log.txt


# Clustering with denoising step
echo Clustering with denoising step
./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance2.biom --output-fasta test/seeds2.fasta --output-compo test/swarms_composition2.tsv \
                --log-file test/log2.txt \
                --denoising

# Clustering with fastidious option
echo Clustering with fastidious option
./clustering.py --input-fasta data/derep.fasta --input-count data/count.tsv \
                --distance 1 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance3.biom --output-fasta test/seeds3.fasta --output-compo test/swarms_composition3.tsv \
                --log-file test/log3.txt \
                --fastidious

# Clustering on its
echo Clustering on its
./clustering.py --input-fasta data/derep_its.fasta --input-count data/count_its.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance4.biom --output-fasta test/seeds4.fasta --output-compo test/swarms_composition4.tsv \
                --log-file test/log4.txt  

# Clustering with denoising step on its
echo Clustering with denoising step on its
./clustering.py --input-fasta data/derep_its.fasta --input-count data/count_its.tsv \
                --distance 2 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance5.biom --output-fasta test/seeds5.fasta --output-compo test/swarms_composition5.tsv \
                --log-file test/log5.txt  \
                --denoising

# Clustering with fastidious opt on its
echo Clustering with fastidious opt on its
./clustering.py --input-fasta data/derep_its.fasta --input-count data/count_its.tsv \
                --distance 1 \
                --nb-cpus 1 \
                --output-biom test/swarms_abundance6.biom --output-fasta test/seeds6.fasta --output-compo test/swarms_composition6.tsv \
                --log-file test/log6.txt  \
                --fastidious
