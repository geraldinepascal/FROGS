#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./FPStep1.py \
 -i data/FPSteps.fasta \
 -b data/FPSteps.biom \
 -p sepp \
 -o test/test_FPStep1.tree \
 -e test/test_FPStep1_excluded.txt \
 -s test/test_FPStep1.fasta \
 -m test/test_FPStep1.biom \
 -c test/test_FPStep1_closests_ref.tsv \
 -t test/test_FPStep1.summary \
 -l test/test_FPStep1.log
