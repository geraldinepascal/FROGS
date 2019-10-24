#!/bin/sh
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$FROGS_DIR/app:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

python deseq2_preprocess.py --data data/data.Rdata \
                            --var EnvType \
                            --log-file test/EnvType_deseq2_preprocess.log \
                            --out-Rdata test/dds_EnvType.rdata

python deseq2_preprocess.py --data data/data.Rdata \
                            --var FoodType \
                            --log-file test/FoodType_deseq2_preprocess.log \
                            --out-Rdata test/dds_FoodType.rdata