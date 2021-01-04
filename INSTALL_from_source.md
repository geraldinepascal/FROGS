We present here guidelines to install dependencies from sources.

It has been tested on a Xubuntu 16.04 virtual machine.

## Installation directories

Here we suppose to install dependencies in the same directory as FROGS.

```bash
version=3.2.0
DIR=`pwd`
BIN_DIR=$DIR/bin
mkdir -p $BIN_DIR
FROGS_libexec=$DIR/FROGS-$version/libexec
FROGS_lib=$DIR/FROGS-$version/lib
FROGS_test=$DIR/FROGS-$version/test
```

## Download FROGS

**download**

```bash
cd $DIR
wget https://github.com/geraldinepascal/FROGS/archive/v$version.tar.gz
tar -xvzf v$version.tar.gz
```

## Installing python and perl dependencies

check python 3 and perl 5 installed (should be done by default on this system )
```bash
python --version
perl --version
```
else 
```bash
sudo apt-get install python3 perl
```

check python package scipy
```python
python
import scipy
# Ctrl + D to quit
```
else
```bash
sudo apt-get install python3-scipy
```



## 1) vsearch 2.15.1, for FROGS Preprocess and FROGS Remove_chimera

**require** :  autoconf, zlib and bzip2 libraries

```bash
sudo apt-get install autoconf libz-dev libbz2-dev
```

**installation**

```bash
cd $BIN_DIR
wget https://github.com/torognes/vsearch/archive/v2.15.1.tar.gz
tar xzf v2.15.1.tar.gz
cd vsearch-2.15.1
./autogen.sh
./configure
make
# test installation
./bin/vsearch -version
# add to FROGS
ln -s $BIN_DIR/vsearch-2.15.1/bin/vsearch $FROGS_libexec/.
```

## 2) FLASH 1.2.11 (optional), for FROGS Preprocess

**installation**

```bash
cd $BIN_DIR
wget https://vorboss.dl.sourceforge.net/project/flashpage/FLASH-1.2.11.tar.gz
tar xvzf FLASH-1.2.11.tar.gz
cd FLASH-1.2.11
make
# check installation
./flash -version
# add to FROGS
ln -s $BIN_DIR/FLASH-1.2.11/flash $FROGS_libexec/.
```

## 3) Pear 0.9.10 (if licence agreement), for FROGS Preprocess

ask for download link and follow installation instructions

## 4) cutadpat 3.1, for FROGS Preprocess

**require** :  pip3
```bash
sudo apt-get install python3-pip
```

**installation**
```bash
cd $BIN_DIR
mkdir cutadapt-3.1
cd cutadapt-3.1
# solution 1 : precise directory (not recommended if you want to use Galaxy)
  sudo pip3 install --prefix=$BIN_DIR/cutadapt-3.1 cutadapt==3.1
  # add cutadapt python library to your PYTHONPATH
  echo export PYTHONPATH="$BIN_DIR/cutadapt-3.1/lib/python3.??/site-packages:\$PYTHONPATH" >> ~/.bashrc
  # check installation
  ./bin/cutadapt --version
  # add to FROGS
  ln -s $BIN_DIR/cutadapt-3.1/bin/cutadapt $FROGS_libexec/.

# solution 2 let pip3 install cutadapt (binary will be available in your PATH)
  #   in your home directory ~/.local/bin
      pip3 install cutadapt==3.1
  #   using sudo in /usr/local/bin
      sudo pip3 install cutadapt==3.1
  # add to FROGS
  link=`which cutadapt`
  ln -s $link $FROGS_libexec/.
```
## 5) swarm 3.0.0, for FROGS Clustering

**installation**
```bash
cd $BIN_DIR
wget https://github.com/torognes/swarm/archive/v3.0.0.tar.gz
tar -vxzf v3.0.0.tar.gz
cd swarm-3.0.0/src
make
cd ../
# check installation
./bin/swarm --version
# add to FROGS
ln -s $BIN_DIR/swarm-3.0.0/bin/swarm $FROGS_libexec/.
```

## 6) ITSx 1.1.2, , for FROGS ITSx

**require** : HMMER >= 3.0

```bash
cd $BIN_DIR
wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
tar xvzf hmmer-3.3.2.tar.gz
cd hmmer-3.3.2/
./configure
make
# check installation
./src/hmmpress -h
./src/hmmscan -h
```

**installation**
```bash
cd $BIN_DIR
wget http://microbiology.se/sw/ITSx_1.1.2.tar.gz
tar -xvzf ITSx_1.1.2.tar.gz
cd ITSx_1.1.2/
# check installation
./ITSx -h 
ln -s $BIN_DIR/ITSx_1.1.2/ITSx $FROGS_libexec/.
ln -s $BIN_DIR/ITSx_1.1.2/ITSx_db $FROGS_libexec/.
```
**recompile ITSx hmm files with your HMMER version**

```bash
cd $BIN_DIR/ITSx_1.1.2/
rm ITSx_db/HMMs/*.h3*
for file in ITSx_db/HMMs/*.hmm
do
$BIN_DIR/hmmer-3.3.2/src/hmmpress $file
done
# this will return an error because of empty N.hmm file, not our fault
```

## 7) NCBI Blast+ blastn 2.10.1, for FROGS Affiliation_OTU

**installation**
```bash
cd $BIN_DIR
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz
tar xvzf ncbi-blast-2.10.1+-x64-linux.tar.gz
# check installation
cd ncbi-blast-2.10.1+/bin/
./blastn -version
# add to FROGS
ln -s $BIN_DIR/ncbi-blast-2.10.1+/bin/blastn $FROGS_libexec/.
```

## 8) RDPClassifier 2.0.3, for FROGS Affiliation_OTU

**require** : git ant and java se jdk
```bash
sudo apt-get install git ant openjdk-8-jdk
```

**installation**
```bash
cd $BIN_DIR
git clone https://github.com/rdpstaff/RDPTools.git
cd RDPTools
git checkout 2.0.3
git submodule init
git submodule update
make
# add to FROGS
ln -s $BIN_DIR/RDPTools/classifier.jar $FROGS_libexec/.
```

## 9) Needlall 6.6.0.0, for FROGS Affiliation_OTU

**require** : pdf and png support
```bash
sudo apt-get install libhpdf-dev libpng-dev
```

**installation**

```bash
cd $BIN_DIR
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar xvf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0 
./configure --prefix=$BIN_DIR/EMBOSS-6.6.0
make 
make install
# check installation
emboss/needleall -h 
# add to FROGS
ln -s $BIN_DIR/EMBOSS-6.6.0/emboss/needleall $FROGS_libexec/.
```

## 10) MAFFT 7.475, for FROGS Tree
**installation**
```bash
cd $BIN_DIR
wget https://mafft.cbrc.jp/alignment/software/mafft-7.475-with-extensions-src.tgz
tar -xvzf mafft-7.475-with-extensions-src.tgz
cd mafft-7.475-with-extensions
```
Edit core/Makefile
* change `PREFIX = /usr/local` with your MAFFT directory (like `# PREFIX = /home/frogs/bin/mafft-7.475-with-extensions/`
```bash
cd core
make clean
make
make install
# check installation
$BIN_DIR/mafft-7.475-with-extensions/bin/mafft -h
# add to FROGS
ln -s $BIN_DIR/mafft-7.475-with-extensions/scripts/mafft $FROGS_libexec/.
```

## 11) FastTree 2.1.11, for FROGS Tree

**installation**
```bash
cd $BIN_DIR
mkdir Fasttree
cd Fasttree
wget http://www.microbesonline.org/fasttree/FastTree # this may change the version!
chmod 777 FastTree
# check install
./FastTree -h
# add to FROGS
ln -s $BIN_DIR/Fasttree/FastTree $FROGS_libexec/.
```

## 12) R 3.6.3, for all FROGSSTAT Phyloseq tools

**installation**

add repository to `/etc/apt/sources.list` for  your ubuntu version! (see https://pbil.univ-lyon1.fr/CRAN/ , this is a french mirror), write `deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/ `at the end of the sources.list file
```bash
# add secure key 
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt-get update
sudo apt-get install r-base
# check installation
R --version
# add to FROGS
link=`which Rscript`
ln -s $link $FROGS_libexec/.
```
### R package dependencies

**require (outside R) **  : httr which need openssl and curl (for plotly)

```bash
sudo apt-get install libssl-dev libcurl4-openssl-dev libfontconfig1-dev libgit2-dev libcairo2-dev

```

**installation**
inside R (use sudo if you want to share installation, else packages will be installed in ~/R):
`R`

* plotly
```R
install.packages("plotly")
# check installation
library(plotly)
```

* phangorn
```R
install.packages("phangorn", dependencies = TRUE)
# check installation
library(phangorn)
```

* phyloseq (this will take some times)
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
# it should install automatically package dependencies
# check installation
library(phyloseq)
```

* DESeq2

  ```R
  BiocManager::install("DESeq2")
  
  # I had to previously (re)install foreign for Hmisc
  # install.packages("foreign")
  # install.packages("Hmisc")
  
  library(DESeq2)
  ```

* optparse

  ```R
  install.packages("optparse")
  # check installation
  library(optparse)
  ```

* calibrate

  ```R
  install.packages("calibrate")
  # check installation
  library(calibrate)
  ```

* formattable

  ```R
  install.packages("formattable")
  # check installation
  library(formattable)
  ```

### validation 

```
sessionInfo()

other attached packages:
 [1] formattable_0.2.0.1         calibrate_1.7.7            
 [3] MASS_7.3-53                 optparse_1.6.6             
 [5] DESeq2_1.26.0               SummarizedExperiment_1.16.1
 [7] DelayedArray_0.12.3         BiocParallel_1.20.1        
 [9] matrixStats_0.57.0          Biobase_2.46.0             
[11] GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
[13] IRanges_2.20.2              S4Vectors_0.24.4           
[15] BiocGenerics_0.32.0         phyloseq_1.30.0            
[17] phangorn_2.5.5              ape_5.4-1                  
[19] plotly_4.9.2.2              ggplot2_3.3.3  
```
versions may change!

### additionnal R and javascript libraries
```
cd $FROGS_lib/external-lib
wget https://github.com/mahendra-mariadassou/phyloseq-extended/archive/v0.99.tar.gz
tar xvf v0.99.tar.gz
mv phyloseq-extended-0.99/* .
rm -r phyloseq-extended-0.99/
```

## 13) pandoc 2.11.3,  for all FROGSSTAT Phyloseq tools
**installation**
```
cd $BIN_DIR
mkdir pandoc
cd pandoc
wget https://github.com/jgm/pandoc/releases/download/2.11.3.2/pandoc-2.11.3.2-1-amd64.deb
sudo dpkg -i pandoc-2.11.3.2-1-amd64.deb
# add to FROGS
link=`which pandoc`
ln -s $link $FROGS_libexec/.
```

## Test FROGS
To check your installation you can type:
```
cd $FROGS_test
# sh test.sh ../ <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
sh test.sh ../ 1 2 results
# Note: JAVA_MEM must be at least 2 (= 2Gb of RAM).
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of Chaillou et al, ISME 2014, doi:10.1038/ismej.2014.202

This test executes the FROGS tools in command line mode.
Example:

```
[user@computer]$ sh test.sh $DIR/FROGS-$version 2 4 res
Step preprocess : Flash lundi 4 janvier 2021, 14:03:15 (UTC+0100)
Step preprocess : Vsearch lundi 4 janvier 2021, 14:06:56 (UTC+0100)
Step clustering lundi 4 janvier 2021, 14:09:29 (UTC+0100)
Step remove_chimera lundi 4 janvier 2021, 14:11:03 (UTC+0100)
Step otu filters lundi 4 janvier 2021, 14:16:18 (UTC+0100)
Step ITSx lundi 4 janvier 2021, 14:16:35 (UTC+0100)
Step affiliation_OTU lundi 4 janvier 2021, 14:16:35 (UTC+0100)
Many-to-many pairwise alignments of two sequence sets
Step affiliation_filter: masking mode lundi 4 janvier 2021, 14:16:43 (UTC+0100)
Step affiliation_filter: deleted mode lundi 4 janvier 2021, 14:16:44 (UTC+0100)
Step affiliation_postprocess lundi 4 janvier 2021, 14:16:45 (UTC+0100)
Step normalisation lundi 4 janvier 2021, 14:16:46 (UTC+0100)
Step clusters_stat lundi 4 janvier 2021, 14:16:46 (UTC+0100)
Step affiliations_stat lundi 4 janvier 2021, 14:16:47 (UTC+0100)
Step biom_to_tsv lundi 4 janvier 2021, 14:16:48 (UTC+0100)
Step biom_to_stdBiom lundi 4 janvier 2021, 14:16:48 (UTC+0100)
Step tsv_to_biom lundi 4 janvier 2021, 14:16:48 (UTC+0100)
Step tree lundi 4 janvier 2021, 14:16:49 (UTC+0100)
Loading required package: ape
Step phyloseq_import_data lundi 4 janvier 2021, 14:17:12 (UTC+0100)
Step phyloseq_composition lundi 4 janvier 2021, 14:18:09 (UTC+0100)
Step phyloseq_alpha_diversity lundi 4 janvier 2021, 14:19:00 (UTC+0100)
Step phyloseq_beta_diversity lundi 4 janvier 2021, 14:19:31 (UTC+0100)
Step phyloseq_structure lundi 4 janvier 2021, 14:19:47 (UTC+0100)
Step phyloseq_clustering lundi 4 janvier 2021, 14:20:38 (UTC+0100)
Step phyloseq_manova lundi 4 janvier 2021, 14:20:54 (UTC+0100)
Step deseq2_preprocess lundi 4 janvier 2021, 14:21:10 (UTC+0100)
Step deseq2_visualization lundi 4 janvier 2021, 14:21:43 (UTC+0100)
Completed with success
```
