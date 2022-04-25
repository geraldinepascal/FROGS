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
./tsv_to_biom.py -t data/swarm_abundance.tsv \
				 -b test/swarm_abundance.biom \
				 --log-file test/swarm_abundance.log

# Without affiliation without sequences (check error)
#./tsv_to_biom.py -t data/swarm_abundance.tsv \
#				 -b test/swarm_abundance.biom \
#				 -f test/swarm_abundance.fasta \
#				 --log-file test/swarm_abundance.log

# Without affiliation with sequences
./tsv_to_biom.py -t data/swarm_fasta_abundance.tsv \
				 -b test/swarm_fasta_abundance.biom \
				 -f test/swarm_fasta_abundance.fasta \
				 --log-file test/swarm_fasta_abundance.log 

# With affiliation with sequences with multihit
./tsv_to_biom.py -t data/affiliation_modifie.tsv \
                 -m data/affiliation_multihit.tsv \
				 -b test/affiliation.biom \
				 -f test/affiliation.fasta \
				 --log-file test/affiliation.log 

# Real dataset
#./tsv_to_biom.py -t /home/maria/Documents/projets/frogs/dev/tsv2biom/TG1101tableFROGScorrigeNom.tsv \
#                 -m /home/maria/Documents/projets/frogs/dev/tsv2biom/TG1101tableFROGS_original_multihit_modified.tsv \
#				 -b test/TG1101.biom \
#				 -f test/TG1101.fasta 