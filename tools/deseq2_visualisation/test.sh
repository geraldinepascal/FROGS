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

# DESeq2 visualisation with OTU 
echo ""
OUT=test/deseq2_otu
echo $OUT "DESeq2 otu abundances"
mkdir -p $OUT

# test on quantitative variables for OTU abundances
./deseq2_visualisation.py --abundanceData  data/phyloseq_OTU.Rdata -a OTU --dds data/EnvType_DESeq_OTU.Rdata \
                            --var EnvType --mod1 SaumonFume --mod2 DesLardons \
                            --log-file $OUT/deseq2_preprocess_EnvType_OTU.log  \
                            --html $OUT/EnvType_DesLardons_SaumonFume_OTU.nb.html 

# DESeq2 visualisation with OTU 
echo ""
OUT=test/deseq2_func
echo $OUT "DESeq2 function abundances"
mkdir -p $OUT

./deseq2_visualisation.py --abundanceData  data/phyloseq_FUNC.Rdata -a FUNC --dds data/EnvType_DESeq_FUNC.Rdata \
                            --var EnvType --mod1 SaumonFume --mod2 DesLardons \
                            --log-file $OUT/deseq2_preprocess_EnvType_FUNC.log  \
                            --html $OUT/EnvType_DesLardons_SaumonFume_FUNC.nb.html 
