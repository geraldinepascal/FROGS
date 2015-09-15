#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

# Blast consensus tax
./affiliations_stat.py \
  --input-biom data/set500_B_affiliation.biom \
  --output-file test/affiliations_metrics1.html \
  --log-file test/log1.txt \
  --tax-consensus-tag "blast_taxonomy" \
  --identity-tag "perc_identity" \
  --coverage-tag "perc_query_coverage" \
  --multiple-tag "blast_affiliations" \
  --rarefaction-ranks Family Genus Species

# Blast tax select first
./affiliations_stat.py \
  --input-biom data/set500_B_affiliation.biom \
  --output-file test/affiliations_metrics2.html \
  --log-file test/log2.txt \
  --taxonomy-tag "taxonomy" \
  --identity-tag "perc_identity" \
  --coverage-tag "perc_query_coverage" \
  --multiple-tag "blast_affiliations" \
  --rarefaction-ranks Family Genus Species

# RDP
./affiliations_stat.py \
  --input-biom data/set500_B_affiliation.biom \
  --output-file test/affiliations_metrics3.html \
  --log-file test/log3.txt \
  --taxonomy-tag "rdp_taxonomy" \
  --bootstrap-tag "rdp_bootstrap" \
  --rarefaction-ranks Family Genus Species