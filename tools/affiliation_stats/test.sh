#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

# Blast consensus tax
echo "# Blast consensus tax"
./affiliation_stats.py \
  --input-biom data/set500_B_affiliation.biom \
  --output-file test/affiliations_metrics1.html \
  --log-file test/log1.txt \
  --tax-consensus-tag "blast_taxonomy" \
  --identity-tag "perc_identity" \
  --coverage-tag "perc_query_coverage" \
  --multiple-tag "blast_affiliations" \
  --rarefaction-ranks Family Genus Species

# RDP
echo "# RDP"
./affiliation_stats.py \
  --input-biom data/set500_B_affiliation.biom \
  --output-file test/affiliations_metrics3.html \
  --log-file test/log3.txt \
  --taxonomy-tag "rdp_taxonomy" \
  --bootstrap-tag "rdp_bootstrap" \
  --rarefaction-ranks Family Genus Species
