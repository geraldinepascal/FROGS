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

python deseq2_visualization.py --phyloseqData data/phyloseq.Rdata --dds data/EnvType_DESeq_dds.Rdata \
                            --var EnvType --mod1 BoeufHache --mod2 SaumonFume \
                            --log-file test/deseq2_preprocess_EnvType.log \
                            --html test/EnvType_BoeufHache_SaumonFume.html

# variables with 2 modes only
python deseq2_visualization.py --phyloseqData data/phyloseq.Rdata --dds data/FoodType_DESeq_dds.Rdata \
                            --var FoodType \
                            --log-file test/deseq2_preprocess_FoodType.log \
                            --html test/FoodType.html

# test on quantitative variables

# test on variables with counfounding effect