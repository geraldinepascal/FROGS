#!/bin/sh
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH


# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

# without normalisation and default ranks
mkdir -p test/phylo_import_default
python r_import_data.py  \
	-b data/chaillou.biom \
	-s data/sample_data.tsv \
	-t data/tree.nwk \
	-d test/phylo_import_default/phylo_import.Rdata -o test/phylo_import_default/phylo_import.html -l test/phylo_import_default/phylo_import.log

# with normalisation and ranks
mkdir -p test/phylo_import_norm_rank
python r_import_data.py  -n -r Royaume Phylum Classe Ordre Famille Genre Espece \
	-b data/chaillou.biom \
	-s data/sample_data.tsv \
	-t data/tree.nwk \
	-d test/phylo_import_norm_rank/phylo_import.Rdata -o test/phylo_import_norm_rank/phylo_import.html -l test/phylo_import_norm_rank/phylo_import.log