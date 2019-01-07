We present here guidelines to install dependancies from sources.

It has been tested on a Xubuntu 16.04 virtual machine.

## Installation directories

Here we suppose to install dependancies in the same directory as FROGS.

```
version=3.1.0
DIR=`pwd`
BIN_DIR=$DIR/bin
mkdir -p $BIN_DIR
FROGS_libexec=$DIR/FROGS-$version/libexec
FROGS_lib=$DIR/FROGS-$version/lib
FROGS_test=$DIR/FROGS-$version/test
```

## Download FROGS

**download**

```
cd $DIR
wget https://github.com/geraldinepascal/FROGS/archive/v$version.tar.gz
tar -xvzf v$version.tar.gz
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

**require** :  autoconf, zlib and bzip2 libraries

`sudo apt-get install autoconf libz-dev libbz2-dev`

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
# add to FROGS
ln -s $BIN_DIR/vsearch-2.9.1/bin/vsearch $FROGS_libexec/.
```

## 2) FLASH 1.2.11 (optional), for FROGS Preprocess

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

**require** :  pip
`sudo apt-get install python-pip`

**installation**
```
cd $BIN_DIR
mkdir cutadapt-1.18
cd cutadapt-1.18
# solution 1 : precise directory (not recommended if you want to use Galaxy)
  sudo pip install --prefix=$BIN_DIR/cutadapt-1.18 cutadapt==1.18
  # add cutadapt python library to your PYTHONPATH
  echo export PYTHONPATH="$BIN_DIR/cutadapt-1.18/lib/python2.7/site-packages:\$PYTHONPATH" >> ~/.bashrc
  # check installation
  ./bin/cutadapt --version
  # add to FROGS
  ln -s $BIN_DIR/cutadapt-1.18/bin/cutadapt $FROGS_libexec/.

# solution 2 let pip install cutadapt (binary will be available in your PATH)
  #   in your home directory ~/.local/bin
      pip install cutadapt==1.18
  #   using sudo in /usr/local/bin
      sudo pip install cutadapt==1.18
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
./bin/swarm --version
# add to FROGS
ln -s $BIN_DIR/swarm-2.2.2/bin/swarm $FROGS_libexec/.
```

## 6) ITSx 1.0.11, , for FROGS ITSx

**require** : HMMER >= 3.0

```
cd $BIN_DIR
wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
tar xvzf hmmer-3.2.1.tar.gz
cd hmmer-3.2.1/
./configure
make
# check installation
./src/hmmpress -h
./src/hmmscan -h
# add to FROGS
ln -s $BIN_DIR/hmmer-3.2.1/src/hmmpress $FROGS_libexec/.
ln -s $BIN_DIR/hmmer-3.2.1/src/hmmscan $FROGS_libexec/.
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
**recompile ITSx hmm files with your HMMER version**

```
cd $BIN_DIR/ITSx_1.0.11/
rm ITSx_db/HMMs/*.h3*
for file in ITSx_db/HMMs/*.hmm
do
$BIN_DIR/hmmer-3.2.1/src/hmmpress $file
done
# this will return an error because of empty N.hmm file, not our fault
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
ln -s $BIN_DIR/ncbi-blast-2.7.1+/bin/blastn $FROGS_libexec/.
```

## 8) RDPClassifier 2.0.2.1, for FROGS Affiliation_OTU

**require** : ant and java se jdk
`sudo apt-get install ant default-jdk`

**installation**
```
cd $BIN_DIR
git clone https://github.com/rdpstaff/RDPTools.git
cd RDPTools
git checkout 2.0.2
git submodule init
git submodule update
make
# add to FROGS
ln -s $BIN_DIR/RDPTools/classifier.jar $FROGS_libexec/.
```

## 9) taskset, for FROGS Affiliation_OTU

**intallation**
`sudo apt-get install util-linux`

## 10) Needlall 6.6.0.0, for FROGS Affiliation_OTU

**require** : pdf and png support
` sudo apt-get install libhpdf-dev libpng-dev` 

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

## 11) MAFFT 7.407, for FROGS Tree
**installation**
```
cd $BIN_DIR
wget https://mafft.cbrc.jp/alignment/software/mafft-7.407-with-extensions-src.tgz
tar -xvzf mafft-7.407-with-extensions-src.tgz
cd mafft-7.407-with-extensions
```
Edit core/Makefile
* change `PREFIX = /usr/local` with your MAFFT directory (like `# PREFIX = /home/frogs/bin/mafft-7.407-with-extensions/`
```
cd core
make clean
make
make install
# check installation
$BIN_DIR/mafft-7.407-with-extensions/bin/mafft -h
# add to FROGS
ln -s $BIN_DIR/mafft-7.407-with-extensions/scripts/mafft $FROGS_libexec/.
```

## 12) FastTree 2.1.10, for FROGS Tree

**installation**
```
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

## 13) R 3.5.1, for all FROGSSTAT Phyloseq tools

**installation**

add repository to `/etc/apt/sources.list` for  your ubuntu version! (see https://pbil.univ-lyon1.fr/CRAN/ ), write `deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/ `at the end of the sources.list file
```
sudo apt-get update
sudo apt-get install r-base
# check installation
R --version
# add to FROGS
link=`which Rscript`
ln -s $link $FROGS_libexec/.
```
### R package dependancies

**require (outside R) **  : httr which need openssl and curl (for plotly)
`sudo apt-get install libssl-dev libcurl4-openssl-dev`

**installation**
inside R (use sudo if you want to share installation, else packages will be installed in ~/R):
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

* phyloseq (this will take some times)
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
# it should install automatically package dependancies
# check installation
library(phyloseq)

```

### validation 
```
sessionInfo()
other attached packages:
[1] phyloseq_1.24.2      gridExtra_2.3        rmarkdown_1.11      
[4] phangorn_2.4.0       ape_5.2              plotly_4.8.0        
[7] ggplot2_3.1.0        BiocInstaller_1.30.0
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

sh test.sh ~/FROGS <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
# Note: JAVA_MEM must be at least 4 (= 4Gb of RAM).
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of Chaillou et al, ISME 2014, doi:10.1038/ismej.2014.202

This test executes the FROGS tools in command line mode.
Example:

```
[user@computer]$ sh test.sh $DIR/FROGS-$version 2 4 res
Step preprocess : Flash jeudi 27 décembre 2018, 10:26:09 (UTC+0100)
Step preprocess : Vsearch jeudi 27 décembre 2018, 10:29:13 (UTC+0100)
Step clustering jeudi 27 décembre 2018, 10:32:33 (UTC+0100)
Step remove_chimera jeudi 27 décembre 2018, 10:41:35 (UTC+0100)
Step filters jeudi 27 décembre 2018, 10:48:16 (UTC+0100)
Step ITSx jeudi 27 décembre 2018, 10:48:29 (UTC+0100)
Step affiliation_OTU jeudi 27 décembre 2018, 10:48:50 (UTC+0100)
Step affiliation_postprocess jeudi 27 décembre 2018, 10:49:33 (UTC+0100)
Step clusters_stat jeudi 27 décembre 2018, 10:49:33 (UTC+0100)
Step affiliations_stat jeudi 27 décembre 2018, 10:49:35 (UTC+0100)
Step biom_to_tsv jeudi 27 décembre 2018, 10:49:38 (UTC+0100)
Step biom_to_stdBiom jeudi 27 décembre 2018, 10:49:39 (UTC+0100)
Step tsv_to_biom jeudi 27 décembre 2018, 10:49:39 (UTC+0100)
Step tree : mafft jeudi 27 décembre 2018, 10:49:39 (UTC+0100)
Step r_import_data jeudi 27 décembre 2018, 10:50:59 (UTC+0100)
Step r_composition jeudi 27 décembre 2018, 10:52:02 (UTC+0100)
Step r_alpha_diversity jeudi 27 décembre 2018, 10:53:01 (UTC+0100)
Step r_beta_diversity jeudi 27 décembre 2018, 10:53:33 (UTC+0100)
Step r_structure jeudi 27 décembre 2018, 10:53:48 (UTC+0100)
Step r_clustering jeudi 27 décembre 2018, 10:54:45 (UTC+0100)
Step r_manova jeudi 27 décembre 2018, 10:54:59 (UTC+0100)
Completed with success
```
