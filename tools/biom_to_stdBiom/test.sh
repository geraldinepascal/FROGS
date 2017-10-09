#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# BIOM without affiliation
./biom_to_stdBiom.py \
 --input-biom data/abundance.biom \
 --output-biom test/abundance_1.biom \
 --log-file test/log_1.txt

# BIOM with affiliation
./biom_to_stdBiom.py \
 --ref blast_taxonomy \
 --input-biom data/affiliation.biom \
 --output-biom test/abundance_blast.biom \
 --log-file test/log_blast.txt

#  # BIOM with affiliation
./biom_to_stdBiom.py \
 --ref rdp_taxonomy \
 --input-biom data/affiliation.biom \
 --output-biom test/abundance_rdp.biom \
 --log-file test/log_rdp.txt
