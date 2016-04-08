#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./clusters_stat.py --input-biom data/swarm.biom \
                   --output-file test/clusters_metrics.html \
                   --log-file test/log.txt
