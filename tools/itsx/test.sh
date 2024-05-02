#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

echo "ITSx ITS1 trimming"
./itsx.py --input-fasta data/input.fasta --input-biom data/input.biom \
    --region ITS1 \
    --output-biom test/itsx.biom --html test/itsx.html --log-file test/itsx.log --output-fasta test/itsx.fasta --output-removed test/removed.fasta

echo "ITSx ckeck ITS only"
./itsx.py --input-fasta data/input.fasta --input-biom data/input.biom \
    --check-its-only --organism-groups F A \
    --output-biom test/itsx.biom --html test/itsx2.html --log-file test/itsx2.log --output-fasta test/itsx2.fasta --output-removed test/removed2.fasta 
