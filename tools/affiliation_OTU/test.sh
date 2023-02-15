#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
   rm test/*
fi

#Blast
echo "#Blast alignment affiliation methods only"
./affiliation_OTU.py --nb-cpus 4 --reference data/db.fasta \
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/Blast_affiliation.biom \
                     --summary test/Blast_affiliation.html --log-file test/Blast_affiliation.log

#Blast and RDP
echo "#Blast alignment and RDP kmer affiliation methods"
./affiliation_OTU.py --reference data/db.fasta --rdp\
                     --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     --output-biom test/Blast_RDP_affiliation.biom \
                     --summary test/Blast_RDP_summary.html --log-file test/Blast_RDP_aff.log

#Blast/Needle
echo "#Blast and Needleall alignment affiliation methods only"
./affiliation_OTU.py --nb-cpus 4 --reference data/ITS1_dataset/Unite_reduced_ITS.fasta \
                     --input-biom data/ITS1_dataset/filters_modified.biom --input-fasta data/ITS1_dataset/filters_modified.fasta \
                     --output-biom test/Blast_Needle_affiliation.biom \
                     --summary test/Blast_Needle_affiliation.html --log-file test/Blast_Needle_affiliation.log

#Blast/Needle & RDP
echo "#Blast and Needleall alignment and RDP kmer affiliation methods"
./affiliation_OTU.py --nb-cpus 4 --rdp --reference data/ITS1_dataset/Unite_reduced_ITS.fasta \
                     --input-biom data/ITS1_dataset/filters_modified.biom --input-fasta data/ITS1_dataset/filters_modified.fasta \
                     --output-biom test/Blast_Needle_affiliation.biom \
                     --summary test/Blast_Needle_RDP_affiliation.html --log-file test/Blast_Needle_RDP_affiliation.log