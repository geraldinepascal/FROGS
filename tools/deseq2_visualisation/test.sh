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

# test on quantitative variables for OTU abundances
python deseq2_visualisation.py --abundanceData  data/phyloseq_OTU.Rdata -a OTU --dds data/EnvType_DESeq_dds.Rdata \
                            --var EnvType --mod1 SaumonFume --mod2 DesLardons \
                            --log-file test/deseq2_preprocess_EnvType_OTU.log  \
                            --html test/EnvType_DesLardons_SaumonFume.nb.html --debug

# test on function abundances
python deseq2_visualisation.py --abundanceData data/

