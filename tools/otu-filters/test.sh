#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
	rm -r test
fi


# FROGS BIOM after affiliation
echo ""
OUT=test/all_filters
echo $OUT
mkdir -p $OUT


# FROGS BIOM before affiliation
echo ""
OUT=test/frogs_woAffiliations
echo $OUT
mkdir -p $OUT
./otu-filters.py \
--input-biom data/set500B_Remove_chimera.biom --input-fasta data/set500B_Remove_chimera.fasta \
--output-biom $OUT/filtered.biom --output-fasta $OUT/filtered.fasta \
--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
--log-file $OUT/log.txt \
--min-sample-presence 2 \
--min-abundance 4296 \
--contaminant data/phi.fa \
--nb-biggest-otu 5


# External BIOM without affiliation
echo ""
OUT=test/usearch_woAffiliations
echo $OUT
mkdir -p $OUT
./otu-filters.py \
--input-biom data/usearch_woAffi.biom --input-fasta data/usearch_woAffi.fasta \
--output-biom $OUT/filtered.biom --output-fasta $OUT/filtered.fasta \
--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
--log-file $OUT/log.txt \
--min-abundance 0.00005 \
--min-sample-presence 2 \
--contaminant data/multi_conta.fa

