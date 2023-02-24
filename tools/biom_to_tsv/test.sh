#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# Without affiliation without sequences
./biom_to_tsv.py --input-biom data/clustering.biom \
                 --output-tsv test/abundance.tsv \
                 --log-file test/log.txt

# Without affiliation with sequences
./biom_to_tsv.py --input-biom data/clustering.biom --input-fasta data/clustering.fasta \
                 --output-tsv test/abundance2.tsv \
                 --log-file test/log2.txt

# With affiliation without sequences
./biom_to_tsv.py --input-biom data/affiliation.biom \
                 --output-tsv test/abundance3.tsv \
                 --log-file test/log3.txt

# With affiliation with sequences
./biom_to_tsv.py --input-biom data/affiliation.biom --input-fasta data/affiliation.fasta \
                 --output-tsv test/abundance4.tsv \
                 --log-file test/log4.txt

# With affiliation without sequences with multi-affi extraction
./biom_to_tsv.py --input-biom data/affiliation.biom \
                 --output-tsv test/abundance5.tsv --output-multi-affi test/multi_affi5.tsv \
                 --log-file test/log5.txt

# With affiliation with sequences with multi-affi extraction
./biom_to_tsv.py --input-biom data/affiliation.biom --input-fasta data/affiliation.fasta \
                 --output-tsv test/abundance6.tsv --output-multi-affi test/multi_affi6.tsv \
                 --log-file test/log6.txt