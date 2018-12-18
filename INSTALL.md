We present here guidelines to install dependancies from sources.

It has been tested on a Xubuntu 16.04 virtual machine.

## Installation directories

Here we suppose to install dependancies in the same directory as FROGS.

```
DIR=`pwd`
BIN_DIR=$DIR/bin
mkdir -p $BIN_DIR
FROGS_libexec=$DIR/FROGS/libexec
FROGS_lib=$DIR/FROGS/lib
FROGS_test=$DIR/FROGS/test
```

## Download FROGS

**require**:  git

`sudo apt install git`

**download**

```
git clone https://github.com/geraldinepascal/FROGS.git
```

## Installing python and perl dependancies

check python 2.7 and perl 5 installed (should be done by default on this system )
```
python --version
perl --version
```
else 
```
sudo apt-get install python perl
```

check python package scipy
```
python
import scipy
# Ctrl + D to quit
```
else
```
sudo apt-get install python-scipy
```



## 1) vsearch 2.9.1, for FROGS Preprocess and FROGS Remove_chimera

**require** :  autoconf
`sudo apt-get install autoconf`

**installation**

```
cd $BIN_DIR
wget https://github.com/torognes/vsearch/archive/v2.9.1.tar.gz
tar xzf v2.9.1.tar.gz
cd vsearch-2.9.1
./autogen.sh
./configure
make
# test installation
./bin/vsearch -version
# add vsearch in FROGS/libexec directory
ln -s $BIN_DIR/vsearch-2.9.1/bin/vsearch $FROGS_libexec/.
```

## 2) FLASH 1.2.11 (optional), for FROGS Preprocess

**require** : zlib
`sudo apt-get install libz-dev`

**installation**

```
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

## 4) cutadpat 1.18, for FROGS Preprocess

**require** :  cython
`sudo apt-get install cython`

**installation**
```
cd $BIN_DIR
wget https://github.com/marcelm/cutadapt/archive/v1.18.tar.gz
tar xvf v1.18.tar.gz
cd cutadapt-1.18
sudo python setup.py install
# check installation
cutadapt --version
# add to FROGS
link=`which cutadapt`
ln -s $link $FROGS_libexec/.
```
## 5) swarm 2.2.2, for FROGS Clustering

**installation**
```
cd $BIN_DIR
wget https://github.com/torognes/swarm/archive/v2.2.2.tar.gz
tar -vxzf v2.2.2.tar.gz
cd swarm-2.2.2/src
make
cd ../
# check installation
./bin/swarm -version
# add to FROGS
ln -s $BIN_DIR/swarm-2.2.2/bin/swarm $FROGS_libexec/.
```

## 6) ITSx 1.0.11, , for FROGS ITSx

**require** HMMER3
```
cd $BIN_DIR
wget http://eddylab.org/software/hmmer/hmmer-3.0.tar.gz
tar xvzf hmmer-3.0.tar.gz
cd hmmer-3.0/
./configure
make
# check installation
./src/hmmpress -h
./src/hmmscan -h
# add to FROGS
ln -s $BIN_DIR/hmmer-3.0/src/hmmpress $FROGS_libexec/.
ln -s $BIN_DIR/hmmer-3.0/src/hmmscan $FROGS_libexec/.
```

**installation**
```
cd $BIN_DIR
wget http://microbiology.se/sw/ITSx_1.0.11.tar.gz
tar -xvzf ITSx_1.0.11.tar.gz
cd ITSx_1.0.11/
# check installation
./ITSx -h 
ln -s $BIN_DIR/ITSx_1.0.11/ITSx $FROGS_libexec/.
ln -s $BIN_DIR/ITSx_1.0.11/ITSx_db $FROGS_libexec/.
```
**recompile ISx hmm files with your HMMER version
```
cd $BIN_DIR/ITSx_1.0.11/
rm ITSx_db/HMMs/*.h3*
for file in ITSx_db/HMMs/*.hmm
do
$BIN_DIR/hmmer-3.0/src/hmmpress $file
done
```

## 7) NCBI Blast+ blastn 2.7.1, for FROGS Affiliation_OTU

**installation**
```
cd $BIN_DIR
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz
tar xvzf ncbi-blast-2.7.1+-x64-linux.tar.gz
# check installation
cd ncbi-blast-2.7.1+/bin/
./blastn -version
# add to FROGS
ln -s $BIN_DIR/ncbi-blast-2.2.30+/bin/blastn $FROGS_libexec/.
```

## 8) RDPClassifier 2.0.2.1, for FROGS Affiliation_OTU

**require** : ant java se jdk
`sudo apt-get install ant default-jdk`

**installation**
```
cd $BIN_DIR
git clone https://github.com/rdpstaff/RDPTools.git
cd RDPTools
git submodule init
git submodule update
make
# add to FROGS
ln -s $BIN_DIR/RDPTools/classifier.jar $FROGS_libexec/.
```

## 9) taskset, for FROGS Affiliation_OTU

** intallation**
`sudo apt-get install util-linux`

## 10) Needlall 6.6.0.0, for FROGS Affiliation_OTU

**require** : pdf, png support
` sudo apt-get install libhpdf-dev libpng-dev

**installation**
```
cd $BIN_DIR
wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
tar xvf EMBOSS-6.6.0.tar.gz
cd EMBOSS-6.6.0 
./configure --prefix=/home/frogs/bin/EMBOSS-6.6.0
make 
make install
# check installation
emboss/needleall -h 
# add to FROGS
ln -s $BIN_DIR/EMBOSS-6.6.0/emboss/needleall $FROGS_libexec/.
```

## 11) MAfft 7.407, , for FROGS Tree
**installation**
```
cd $BIN_DIR
wget https://mafft.cbrc.jp/alignment/software/mafft-7.407-with-extensions-src.tgz
tar -xvzf mafft-7.407-with-extensions-src.tgz
cd mafft-7.407-with-extensions/
cd core
```
Edit Makefile
* change `PREFIX = /usr/local` with your $DIR (like `# PREFIX = /home/frogs/bin/mafft-7.407-with-extensions/`
* change `BINDIR = $(PREFIX)/bin` with your $DIR (like `# BINDIR = /home/frogs/bin/mafft-7.407-with-extensions/bin`

```
make clean
make
# check installation
$BIN_DIR/mafft-7.310-with-extensions/bin/mafft -h
# add to FROGS
ln -s $BIN_DIR/mafft-7.310-with-extensions/bin/mafft $FROGS_libexec/.
```

## 12) FastTree 2.1.10, for FROGS Tree

**installation**
```
$BIN_DIR
mkdir Fasttree
cd Fasttree
wget http://www.microbesonline.org/fasttree/FastTree # this may change the version!
chmod 777 FastTree
# check install
FastTree -h
# add to FROGS
ln -s $BIN_DIR/Fasttree/FastTree $FROGS_libexec/.
```

## 13) R 3.5.1, for all FROGSSTAT Phyloseq tools

**installation**
```
# add repository to /etc/apt/sources.list :
sudo echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
sudo apt-get update
sudo apt-get install r-base
# check installation
R --version
# add to FROGS
link=`which Rscript`
ln -s $link FROGS_libexec/.
```
### R package dependancies

**require (outside R) **  : httr which need openssl and curl (for plotly)
`sudo apt-get install libssl libcurl4-openssl-dev`

**installation**
inside R:
`R`

* plotly
```
install.packages("plotly")
# check installation
library(plotly)
```

* phangorn
```
install.packages("phangorn")
# check installation
library(phangorn)
```
* rmarkdown
```
install.packages("rmarkdown")
# check installation
library(rmarkdown)
```

* gridextra
```
install.packages("gridExtra")
# check installation
library(gridExtra)
```

* phyloseq
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
# check installation
library(phyloseq)
```

### validation
```
sessionInfo()
# [1] phyloseq_1.22.3      BiocInstaller_1.28.0 phangorn_2.4.0      
# [4] ape_5.2              rmarkdown_1.10       gridExtra_2.3       
# [7] plotly_4.8.0         ggplot2_3.1.0 
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

## 14) pandoc 1.19.2.1,  for all FROGSSTAT Phyloseq tools
**installation**
```
cd $BIN_DIR
mkdir pandoc
cd pandoc
wget https://github.com/jgm/pandoc/releases/download/1.19.2.1/pandoc-1.19.2.1-1-amd64.deb
sudo dpkg -i pandoc-1.19.2.1-1-amd64.deb
# add to FROGS
ln -s /usr/bin/pandoc $FROGS_libexec/.
```

## Test FROGS
To check your installation you can type:
```
cd $FROGS_test
# Note: JAVA_MEM must be at least 4 (= 4Gb of RAM).

sh test.sh ~/FROGS <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of Chaillou et al, ISME 2014, doi:10.1038/ismej.2014.202

This test executes the FROGS tools in command line mode.
Example:
```
[user@computer:/home/frogs/FROGS/test/]$ sh test.sh ~/FROGS 2 4 res
Step preprocess : Flash mercredi 10 octobre 2018, 14:11:30 (UTC+0200)
Step preprocess : Vsearch mercredi 10 octobre 2018, 14:13:33 (UTC+0200)
Step clustering mercredi 10 octobre 2018, 14:15:36 (UTC+0200)
Step remove_chimera mercredi 10 octobre 2018, 14:18:43 (UTC+0200)
Step filters mercredi 10 octobre 2018, 14:22:36 (UTC+0200)
Step ITSx mercredi 10 octobre 2018, 14:22:42 (UTC+0200)
Step affiliation_OTU mercredi 10 octobre 2018, 14:22:42 (UTC+0200)
Step affiliation_postprocess mercredi 10 octobre 2018, 14:23:08 (UTC+0200)
Step clusters_stat mercredi 10 octobre 2018, 14:23:08 (UTC+0200)
Step affiliations_stat mercredi 10 octobre 2018, 14:23:09 (UTC+0200)
Step biom_to_tsv mercredi 10 octobre 2018, 14:23:12 (UTC+0200)
Step biom_to_stdBiom mercredi 10 octobre 2018, 14:23:12 (UTC+0200)
Step tsv_to_biom mercredi 10 octobre 2018, 14:23:12 (UTC+0200)
Step tree : mafft mercredi 10 octobre 2018, 14:23:26 (UTC+0200)
Step r_import_data mercredi 10 octobre 2018, 14:25:25 (UTC+0200)
Step r_composition mercredi 10 octobre 2018, 14:25:39 (UTC+0200)
Step r_alpha_diversity mercredi 10 octobre 2018, 14:25:53 (UTC+0200)
Step r_beta_diversity mercredi 10 octobre 2018, 14:26:19 (UTC+0200)
Step r_structure mercredi 10 octobre 2018, 14:26:31 (UTC+0200)
Step r_clustering mercredi 10 octobre 2018, 14:26:47 (UTC+0200)
Step r_manova mercredi 10 octobre 2018, 14:26:57 (UTC+0200)
Completed with success
```