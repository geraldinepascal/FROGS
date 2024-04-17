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

# Illumina R1 and R2 vsearch swarm
echo "Illumina R1 and R2 with vsearch and swarm"
./denoising.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software vsearch \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --process swarm \
                         --output-fasta test/derep_illumina_R1R2.fasta --output-biom test/illumina_R1R2.biom --html test/summary_illumina_R1R2.html --log-file test/log_illumina_R1R2.txt

# Illumina R1 and R2 vsearch dada2
echo "Illumina R1 and R2 with vsearch and dada2"
./denoising.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software vsearch \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --process dada2 \
                         --output-fasta test/derep_illumina_dada2_R1R2.fasta --output-biom test/illumina_dada2_R1R2.biom --html test/summary_illumina_dada2_R1R2.html --log-file test/log_illumina_dada2_R1R2.txt

# Illumina R1 and R2 with PEAR swarm
echo "Illumina R1 and R2 with PEAR and swarm"
./denoising.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software pear \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --process swarm \
                         --output-fasta test/derep_illumina_pear_R1R2.fasta --output-biom test/illumina_pear_R1R2.biom --html test/summary_illumina_pear_R1R2.html --log-file test/log_illumina_pear_R1R2.txt

# Illumina R1 and R2 with flash swarm
echo "Illumina R1 and R2 with flash and swarm"
./denoising.py illumina --input-R1 data/sampleA_R1.fastq.gz data/sampleB_R1.fastq.gz --input-R2 data/sampleA_R2.fastq.gz data/sampleB_R2.fastq.gz \
                         --samples-names sample_A sample_B \
                         --R1-size 251 --R2-size 251 \
                         --merge-software flash \
                         --expected-amplicon-size 410 \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --process swarm \
                         --output-fasta test/derep_illumina_flash_R1R2.fasta --output-biom test/illumina_flash_R1R2.biom --html test/summary_illumina_flash_R1R2.html --log-file test/log_illumina_flash_R1R2.txt

# Illumina tar R1 and R2, keep unmerged swarm
echo "Illumina tar R1 and R2, keep unmerged"
tar -zcf test/samples.tar.gz -C data sampleA_R1.fastq.gz sampleA_R2.fastq.gz sampleB_R1.fastq.gz sampleB_R2.fastq.gz
./denoising.py illumina --input-archive test/samples.tar.gz \
                         --R1-size 251 --R2-size 251 --keep-unmerged \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --process swarm \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-fasta test/derep_illumina_TAR_R1R2_keepUnmerged.fasta --output-biom test/count_illumina_TAR_R1R2_keepUnmerged.biom --html test/summary_illumina_TAR_R1R2_keepUnmerged.html --log-file test/log_illumina_TAR_R1R2_keepUnmerged.txt
rm test/samples.tar.gz

# Illumina tar R1 and R2, keep unmerged dada2
echo "Illumina tar R1 and R2, keep unmerged"
tar -zcf test/samples.tar.gz -C data sampleA_R1.fastq.gz sampleA_R2.fastq.gz sampleB_R1.fastq.gz sampleB_R2.fastq.gz
./denoising.py illumina --input-archive test/samples.tar.gz \
                         --R1-size 251 --R2-size 251 --keep-unmerged \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --process dada2 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-fasta test/derep_illumina_dada2_TAR_R1R2_keepUnmerged.fasta --output-biom test/count_illumina_dada2_TAR_R1R2_keepUnmerged.biom --html test/summary_illumina_dada2_TAR_R1R2_keepUnmerged.html --log-file test/log_illumina_dada2_TAR_R1R2_keepUnmerged.txt
rm test/samples.tar.gz

# Illumina contiged
echo "Illumina contiged"
./denoising.py illumina --input-R1 data/sampleA.fastq.gz data/sampleB.fastq.gz --already-contiged \
                         --samples-names sample_A sample_B \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-fasta test/derep_illumina_CONTIGED.fasta --output-biom test/count_illumina_CONTIGED.biom --html test/summary_illumina_CONTIGED.html --log-file test/log_illumina_CONTIGED.txt

# Illumina tar contiged
echo "Illumina tar contiged"
tar -zcf test/contiged_samples.tar.gz -C data sampleA.fastq.gz sampleB.fastq.gz
./denoising.py illumina --input-archive test/contiged_samples.tar.gz --already-contiged \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --five-prim-primer "CCGTCAATTC" --three-prim-primer "CCGCNGCTGCT" \
                         --output-fasta test/derep_illumina_TAR_CONTIGED.fasta --output-biom test/count_illumina_TAR_CONTIGED.tsv --html test/summary_illumina_TAR_CONTIGED.html --log-file test/log_illumina_TAR_CONTIGED.txt
rm test/contiged_samples.tar.gz

# Illumina contiged without primers
echo "Illumina contiged without primers"
cutadapt -g CCGTCAATTC --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 9 -o test/test_sampleA_tmp.fastq.gz data/sampleA.fastq.gz > /dev/null 2>&1
cutadapt -a CCGCNGCTGCT --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 10 -o test/test_sampleA.fastq.gz test/test_sampleA_tmp.fastq.gz > /dev/null 2>&1
cutadapt -g CCGTCAATTC --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 9 -o test/test_sampleB_tmp.fastq.gz data/sampleB.fastq.gz > /dev/null 2>&1
cutadapt -a CCGCNGCTGCT --error-rate 0.1 --discard-untrimmed --match-read-wildcards --minimum-length 200 --overlap 10 -o test/test_sampleB.fastq.gz test/test_sampleB_tmp.fastq.gz > /dev/null 2>&1
./denoising.py illumina --input-R1 test/test_sampleA.fastq.gz test/test_sampleB.fastq.gz \
                         --samples-names sample_A sample_B \
                         --already-contiged --without-primers \
                         --min-amplicon-size 340 --max-amplicon-size 450 \
                         --output-fasta test/derep_illumina_contiged_custom.fasta --output-biom test/count_illumina_contiged_custom.biom --html test/summary_illumina_contiged_custom.html --log-file test/log_illumina_contiged_custom.txt
rm test/test_sampleA.fastq.gz test/test_sampleA_tmp.fastq.gz test/test_sampleB.fastq.gz test/test_sampleB_tmp.fastq.gz

# 454 R1
echo "454 R1"
./denoising.py 454 --input-R1 data/SRR443364_clipped.fastq.gz \
                    --min-amplicon-size 340 --max-amplicon-size 450 \
                    --five-prim-primer "ACGGGAGGCAGCAG" --three-prim-primer "AGGATTAGATACCCTGGTA" \
                    --output-fasta test/derep_SRR443364_454.fasta --output-biom test/count_SRR443364_454.biom --html test/summary_SRR443364_454.html --log-file test/log_SRR443364_454.txt
                    
# Long reads PacBio swarm
echo "Long reads swarm"
./denoising.py longreads --input-archive data/LongReads.tar.gz \
                    --min-amplicon-size 400 --max-amplicon-size 3000 \
                    --five-prim-primer AGRGTTYGATYMTGGCTCAG --three-prim-primer AAGTCGTAACAAGGTARCY \
                    --process swarm \
                    --nb-cpus 4 \
                    --log-file test/denoising_longreads_swarm.log --output-fasta test/denoising_longreads_swarm.fasta --output-biom test/denoising_longreads_swarm.biom --html test/denoising_longreads_swarm.html
			
# Long reads PacBio dada2
echo "Long reads dada2"
./denoising.py longreads --input-archive data/LongReads.tar.gz \
                    --min-amplicon-size 400 --max-amplicon-size 3000 \
                    --five-prim-primer AGRGTTYGATYMTGGCTCAG --three-prim-primer AAGTCGTAACAAGGTARCY \
                    --process dada2 \
                    --nb-cpus 4 \
                    --log-file test/denoising_longreads_dada2.log --output-fasta test/denoising_longreads_dada2.fasta --output-biom test/denoising_longreads_dada2.biom --html test/denoising_longreads_dada2.html
