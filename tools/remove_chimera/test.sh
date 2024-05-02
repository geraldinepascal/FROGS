#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# Input BIOM
./remove_chimera.py \
                    --input-fasta data/swarm_seed_sequences.fasta --input-biom data/swarm_abundance.biom \
                    --output-fasta test/non_chimera.fasta --output-biom test/abundance.biom --html test/summary.html --log-file test/log.txt \
                    --nb-cpus 2
