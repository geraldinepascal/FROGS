#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$FROGS_DIR/app:$PATH
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
OUT=test/metrics-filter
echo $OUT "delete mode"
mkdir -p $OUT
./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_OTU_deleted.biom --output-fasta $OUT/filtered_OTU_deleted.fasta \
--summary $OUT/summary_OTU_deleted.html --impacted $OUT/impacted_OTU_deleted.tsv --impacted-multihit $OUT/impacted_OTU_deleted_multihit.tsv \
--log-file $OUT/delete_log.txt \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--delete

echo ""
echo $OUT "masking mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_OTU_masked.biom --output-fasta $OUT/filtered_OTU_masked.fasta \
--summary $OUT/summary_OTU_masked.html --impacted $OUT/impacted_OTU_masked.tsv --impacted-multihit $OUT/impacted_OTU_masked_multihit.tsv \
--log-file $OUT/mask_log.txt \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--mask

echo ""
OUT=test/taxon-ignored-filter
echo $OUT "delete mode"
mkdir -p $OUT

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_OTU_deleted.biom --output-fasta $OUT/filtered_OTU_deleted.fasta \
--summary $OUT/summary_OTU_deleted.html --impacted $OUT/impacted_OTU_deleted.tsv --impacted-multihit $OUT/impacted_OTU_deleted_multihit.tsv \
--log-file $OUT/delete_log.txt \
--taxon-ignored "Methylovulum miyakonense" "subsp." "unknown species" \
--delete

echo ""
echo $OUT "masking mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_OTU_masked.biom --output-fasta $OUT/filtered_OTU_masked.fasta \
--summary $OUT/summary_OTU_masked.html --impacted $OUT/impacted_OTU_masked.tsv --impacted-multihit $OUT/impacted_OTU_masked_multihit.tsv \
--log-file $OUT/mask_log.txt \
--taxon-ignored "Methylovulum miyakonense" "subsp." "unknown species" \
--mask

OUT=test/all-filter
mkdir -p $OUT
echo ""
echo $OUT "masking mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_OTU_masked.biom --output-fasta $OUT/filtered_OTU_masked.fasta \
--summary $OUT/summary_OTU_masked.html --impacted $OUT/impacted_OTU_masked.tsv --impacted-multihit $OUT/impacted_OTU_masked_multihit.tsv \
--log-file $OUT/mask_log.txt \
--mask \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species --min-rdp-bootstrap Species:0.7303 \
--max-blast-evalue 0 --min-blast-identity 1.0 --min-blast-coverage 1.0 \
--taxon-ignored "Methylovulum miyakonense" "subsp." "unknown species" --debug
