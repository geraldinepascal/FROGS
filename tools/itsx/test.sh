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

echo "ITSx"
./itsx.py --input-fasta data/input.fasta --input-biom data/input.biom --region ITS1 --out-abundance test/itsx.biom --summary test/itsx.html --log-file test/itsx.log --out-fasta test/itsx.fasta
