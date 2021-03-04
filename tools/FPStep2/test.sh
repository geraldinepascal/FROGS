#!/bin/sh

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

#python2.7 gene_placement.py -i 16S -t out5.tree -o marker.tsv.gz -n
# Dans mon ordi
python3 ../FPStep2.py -c 16S -i EC -t sortti_etape1.tre -o sorti_test_EC__etape2.tsv

