#!/bin/bash
export PATH=../../bin:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH
./clusters_stat.py --input-biom data/swarm.biom \
                   --output-file test/clusters_metrics.html \
                   --log-file test/log.txt