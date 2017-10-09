#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

./clusters_stat.py --input-biom data/swarm.biom \
                   --output-file test/clusters_metrics.html \
                   --log-file test/log.txt
