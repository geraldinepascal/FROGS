#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

#Blast
echo "#Blast"
./affiliation_OTU.py --reference data/db.fasta \
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/Blast_affiliation.biom \
                     --summary test/Blast_summary.html --log-file test/Blast_aff.log
#Blast and RDP
echo "#Blast and RDP"
./affiliation_OTU.py --reference data/db.fasta --rdp\
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/Blast_RDP_affiliation.biom \
                     --summary test/Blast_RDP_summary.html --log-file test/Blast_RDP_aff.log
