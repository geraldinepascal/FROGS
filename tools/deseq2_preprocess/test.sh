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

# test on variable with multiple mode
echo "# test on variable with multiple mode"
./deseq2_preprocess.py --data data/data.Rdata \
                            --var EnvType \
                            --log-file test/EnvType_deseq2_preprocess.log \
                            --out-Rdata test/dds_EnvType.rdata

# variables with 2 modes only
echo "# variables with 2 modes only"
./deseq2_preprocess.py --data data/data.Rdata \
                            --var FoodType \
                            --log-file test/FoodType_deseq2_preprocess.log \
                            --out-Rdata test/dds_FoodType.rdata
