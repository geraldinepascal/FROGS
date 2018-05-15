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
./affiliation_OTU.py --nb-cpus 4 --reference data/db.fasta \
                     --input-biom data/filter_modified.biom --input-fasta data/filter_modified.fasta \
                     --output-biom test/Needle_affiliation.biom \
                     --summary test/Needle_summary.html --log-file test/Needle_aff.log
#Blast and RDP
#~ echo "#Blast and RDP"
#~ ./affiliation_OTU.py --reference data/db.fasta --rdp\
                     #~ --input-biom data/swarm.biom --input-fasta data/swarm.fasta \
                     #~ --output-biom test/Blast_RDP_affiliation.biom \
                     #~ --summary test/Blast_RDP_summary.html --log-file test/Blast_RDP_aff.log
