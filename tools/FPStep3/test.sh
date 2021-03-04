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
# Dans mon ordi (cd DonneesPourTest)
python3 ../FPStep3.py -i test_16S.biom1 -m marker.tsv -f out_reaction.tsv --function_abund predmeta.tsv --seqtab seqtab_norm.tsv --weighted weigth

