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
--output-biom $OUT/filtered_ASV_deleted.biom --output-fasta $OUT/filtered_ASV_deleted.fasta \
--summary $OUT/summary_ASV_deleted.html --impacted $OUT/impacted_ASV_deleted.tsv --impacted-multihit $OUT/impacted_ASV_deleted_multihit.tsv \
--log-file $OUT/delete_log.txt \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--delete

echo ""
OUT=test/metrics-filter
echo $OUT "masking mode"
mkdir -p $OUT

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_ASV_masked.biom --output-fasta $OUT/filtered_ASV_masked.fasta \
--summary $OUT/summary_ASV_masked.html --impacted $OUT/impacted_ASV_masked.tsv --impacted-multihit $OUT/impacted_ASV_masked_multihit.tsv \
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
--output-biom $OUT/filtered_ASV_deleted.biom --output-fasta $OUT/filtered_ASV_deleted.fasta \
--summary $OUT/summary_ASV_deleted.html --impacted $OUT/impacted_ASV_deleted.tsv --impacted-multihit $OUT/impacted_ASV_deleted_multihit.tsv \
--log-file $OUT/delete_log.txt \
--ignore-blast-taxa "Methylovulum miyakonense" "subsp." "unknown species" \
--delete

echo ""
echo $OUT "masking mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_ASV_masked.biom --output-fasta $OUT/filtered_ASV_masked.fasta \
--summary $OUT/summary_ASV_masked.html --impacted $OUT/impacted_ASV_masked.tsv --impacted-multihit $OUT/impacted_ASV_masked_multihit.tsv \
--log-file $OUT/mask_log.txt \
--ignore-blast-taxa "Methylovulum miyakonense" "subsp." "unknown species" \
--mask

OUT=test/all-filter
mkdir -p $OUT
echo ""
echo $OUT "deleting mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_ASV_deleted.biom --output-fasta $OUT/filtered_ASV_deleted.fasta \
--summary $OUT/summary_ASV_deleted.html --impacted $OUT/impacted_ASV_deleted.tsv --impacted-multihit $OUT/impacted_ASV_deleted_multihit.tsv \
--log-file $OUT/del_log.txt \
--delete \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--ignore-blast-taxa "Methylovulum miyakonense" "subsp." "unknown species"

echo ""
echo $OUT "masking mode"

./affiliation_filters.py \
--input-biom data/fake_affiliation.biom --input-fasta data/fake_affiliation.fasta \
--output-biom $OUT/filtered_ASV_masked.biom --output-fasta $OUT/filtered_ASV_masked.fasta \
--summary $OUT/summary_ASV_masked.html --impacted $OUT/impacted_ASV_masked.tsv --impacted-multihit $OUT/impacted_ASV_masked_multihit.tsv \
--log-file $OUT/mask_log.txt \
--mask \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 402 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 0 \
--ignore-blast-taxa "Methylovulum miyakonense" "subsp." "unknown species"

OUT=test/without_affiliation
mkdir -p $OUT
echo ""
echo $OUT "delete mode"

./affiliation_filters.py \
--input-biom data/without_affiliation.biom --input-fasta data/without_affiliation.fasta \
--output-biom $OUT/filtered_ASV_masked.biom --output-fasta $OUT/filtered_ASV_masked.fasta \
--summary $OUT/summary_ASV_without_affi_deleted.html --impacted $OUT/summary_ASV_without_affi_deleted.tsv --impacted-multihit $OUT/summary_ASV_without_affi_deleted.tsv \
--log-file $OUT/summary_ASV_without_affi_deleted.txt \
--delete \
--min-blast-length 40 \
--min-blast-identity 0.7 \
--min-blast-coverage 0.2 \
--keep-blast-taxa "Saccharomycetales"