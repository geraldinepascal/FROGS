#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi


# FROGS BIOM after affiliation
echo ""
OUT=test/all_filters
echo $OUT
mkdir -p $OUT
python2.7 filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered.biom --output-fasta $OUT/filtered.fasta \
--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
--log-file $OUT/log.txt \
--min-sample-presence 2 \
--min-abundance 4296 \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--contaminant data/phi.fa \
--nb-biggest-otu 5


# FROGS BIOM before affiliation
echo ""
OUT=test/frogs_woAffiliations
echo $OUT
mkdir -p $OUT
python2.7 filters.py \
--input-biom data/set500B_Remove_chimera.biom --input-fasta data/set500B_Remove_chimera.fasta \
--output-biom $OUT/filtered.biom --output-fasta $OUT/filtered.fasta \
--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
--log-file $OUT/log.txt \
--min-abundance 0.00005 \
--min-sample-presence 2


# External BIOM without affiliation
echo ""
OUT=test/usearch_woAffiliations
echo $OUT
mkdir -p $OUT
python2.7 filters.py \
--input-biom data/usearch_woAffi.biom --input-fasta data/usearch_woAffi.fasta \
--output-biom $OUT/filtered.biom --output-fasta $OUT/filtered.fasta \
--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
--log-file $OUT/log.txt \
--min-abundance 0.00005 \
--min-sample-presence 2


# Large FROGS BIOM
#echo ""
#OUT=test/large_dataset
#echo $OUT
#mkdir -p $OUT
#python2.7 filters.py --input-biom data/large.biom --input-fasta data/large.fasta \
#--output-fasta $OUT/fasta_result.fasta --output-biom $OUT/biom_output \
#--summary $OUT/summary.html --excluded $OUT/excluded.tsv \
#--log-file $OUT/log.txt \
#--min-sample-presence 2 \
#--min-abundance 200
