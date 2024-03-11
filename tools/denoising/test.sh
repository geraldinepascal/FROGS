#!/bin/bash
FROGS_DIR=`dirname $(dirname $(pwd))`
export PATH=$FROGS_DIR/libexec:$PATH
export PYTHONPATH=$FROGS_DIR/lib:$PYTHONPATH

# Create output folder
if [ ! -d "test" ]
then
    mkdir test
else
    rm -r test/*
fi

# Illumina R1 and R2
echo "Illumina R1 and R2 with vsearch"
./preprocess.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software vsearch \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_R1R2.fasta --output-count test/count_illumina_R1R2.tsv --summary test/summary_illumina_R1R2.html --log-file test/log_illumina_R1R2.txt

# Illumina R1 and R2 with PEAR
echo "Illumina R1 and R2 with PEAR"
./preprocess.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software pear \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_R1R2_pear.fasta --output-count test/count_illumina_R1R2_pear.tsv --summary test/summary_illumina_R1R2_pear.html --log-file test/log_illumina_R1R2_pear.txt

# Illumina R1 and R2 with flash
echo "Illumina R1 and R2 with flash"
./preprocess.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --expected-amplicon-size 410 \
                         --merge-software flash \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_R1R2_flash.fasta --output-count test/count_illumina_R1R2_flash.tsv --summary test/summary_illumina_R1R2_flash.html --log-file test/log_illumina_R1R2_flash.txt

# Illumina tar R1 and R2, keep unmerged
echo "Illumina tar R1 and R2, keep unmerged"
tar -zcf test/samples.tar.gz -C data sampleA_R1.fastq.gz sampleA_R2.fastq.gz sampleB_R1.fastq.gz sampleB_R2.fastq.gz
./preprocess.py illumina --input-archive test/samples.tar.gz \
                         --R1-size 251 --R2-size 251 --keep-unmerged \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_TAR_R1R2_keepUnmerged.fasta --output-count test/count_illumina_TAR_R1R2_keepUnmerged.tsv --summary test/summary_illumina_TAR_R1R2_keepUnmerged.html --log-file test/log_illumina_TAR_R1R2_keepUnmerged.txt
rm test/samples.tar.gz

# Illumina contiged
echo "Illumina contiged"
./preprocess.py illumina --input-R1 data/sampleA.fastq.gz data/sampleB.fastq.gz --already-contiged \
                         --samples-names sample_A sample_B \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_CONTIGED.fasta --output-count test/count_illumina_CONTIGED.tsv --summary test/summary_illumina_CONTIGED.html --log-file test/log_illumina_CONTIGED.txt

# Illumina tar contiged
echo "Illumina tar contiged"
tar -zcf test/contiged_samples.tar.gz -C data sampleA.fastq.gz sampleB.fastq.gz
./preprocess.py illumina --input-archive test/contiged_samples.tar.gz --already-contiged \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-dereplicated test/derep_illumina_TAR_CONTIGED.fasta --output-count test/count_illumina_TAR_CONTIGED.tsv --summary test/summary_illumina_TAR_CONTIGED.html --log-file test/log_illumina_TAR_CONTIGED.txt
rm test/contiged_samples.tar.gz

# Illumina contiged without primers
echo "Illumina contiged without primers"
cutadapt -g CCGTCAATTC --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 9 -o test/test_sampleA_tmp.fastq.gz data/sampleA.fastq.gz > /dev/null 2>&1
cutadapt -a CCGCNGCTGCT --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 10 -o test/test_sampleA.fastq.gz test/test_sampleA_tmp.fastq.gz > /dev/null 2>&1
cutadapt -g CCGTCAATTC --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 9 -o test/test_sampleB_tmp.fastq.gz data/sampleB.fastq.gz > /dev/null 2>&1
cutadapt -a CCGCNGCTGCT --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 10 -o test/test_sampleB.fastq.gz test/test_sampleB_tmp.fastq.gz > /dev/null 2>&1
./preprocess.py illumina --input-R1 test/test_sampleA.fastq.gz test/test_sampleB.fastq.gz \
                         --samples-names sample_A sample_B \
                         --already-contiged --without-primers \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --output-dereplicated test/derep_illumina_contiged_custom.fasta --output-count test/count_illumina_contiged_custom.tsv --summary test/summary_illumina_contiged_custom.html --log-file test/log_illumina_contiged_custom.txt
rm test/test_sampleA.fastq.gz test/test_sampleA_tmp.fastq.gz test/test_sampleB.fastq.gz test/test_sampleB_tmp.fastq.gz

# 454 R1
echo "454 R1"
./preprocess.py 454 --input-R1 data/SRR443364_clipped.fastq.gz \
                    --min-amplicon-size 340 --max-amplicon-size 450 \
                    --five-prim-primer "ACGGGAGGCAGCAG" --three-prim-primer "AGGATTAGATACCCTGGTA" \
                    --output-dereplicated test/derep_SRR443364_454.fasta --output-count test/count_SRR443364_454.tsv --summary test/summary_SRR443364_454.html --log-file test/log_SRR443364_454.txt
