#!/bin/bash
export PATH=../../libexec:$PATH
export PYTHONPATH=../../bin:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
fi

OUT=test/test_pe/both
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py \
  --input-R1 data/test2_R1.fq --input-R2 data/test2_R2.fq --input-barcode data/barcode2.txt \
  --mismatches 1 --end both \
  --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt 
echo ""

OUT=test/test_pe/forward
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/test2_R1.fq --input-R2 data/test2_R2.fq --input-barcode data/barcode2_forward.txt --mismatches 1 --end bol\
                 --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt 
echo ""

OUT=test/test_pe/reverse
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/test2_R1.fq --input-R2 data/test2_R2.fq --input-barcode data/barcode2_reverse.txt --mismatches 1 --end eol\
                 --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt                  
echo ""
 
OUT=test/test_se/both
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/test2_R1.fq --input-barcode data/barcode2.txt --mismatches 1 --end both\
  --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt 
echo ""

OUT=test/test_se/forward
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/test2_R1.fq --input-barcode data/barcode2_forward.txt --mismatches 1 --end bol\
  --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt 
echo ""

OUT=test/test_se/reverse
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/test2_R1.fq --input-barcode data/barcode2_reverse.txt --mismatches 1 --end eol\
                 --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt                  
echo ""

OUT=test/formation_se/both
mkdir -p $OUT
echo $OUT
rm -r $OUT/*
./demultiplex.py --input-R1 data/formation.fastq --input-barcode data/formation.barcode.txt --mismatches 1 --end both\
                 --output-demultiplexed $OUT/demultiplexed.tar.gz --output-excluded $OUT/undemultiplexed.tar.gz --log-file $OUT/log.txt --summary $OUT/summary.txt 
echo ""
