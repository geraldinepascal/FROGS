<p align="center">
 <a href="http://frogs.toulouse.inra.fr/">
  <img src="FROGS_logo.png" align="center" width="20%" style="display: block; margin: auto;"/>
 </a>
</p>

Visit our web site : http://frogs.toulouse.inrae.fr/


[![Release](https://img.shields.io/badge/release-3.2.0-blue.svg)![Date](https://img.shields.io/badge/date-May%202020-red.svg)](https://github.com/geraldinepascal/FROGS-wrappers/releases) [<img src="https://www.podcastscience.fm/wp-content/uploads/2017/12/deezer.png" width="5%" style="display: block; margin: auto;"/>](https://www.deezer.com/fr/playlist/5233843102?utm_source=deezer&utm_content=playlist-5233843102&utm_term=18632989_1545296531&utm_medium=web)



# Description

FROGS is a CLI workflow designed to produce an OTU count matrix from high depth sequencing amplicon data.

FROGS-wrappers allow to add FROGS on a Galaxy instance. (see https://github.com/geraldinepascal/FROGS-wrappers)

This workflow is focused on:

- User-friendliness with lots of rich graphic outputs and the integration in Galaxy thanks to FROGS-wrappers
- Accuracy with a clustering without global similarity threshold, the management of multi-affiliations and management of separated PCRs in the chimera removal step
- Speed with fast algorithms and an easy to use parallelisation
- Scalability with algorithms designed to support the data growth



# Table of content

* [Convenient input data](#convenient-input-data)
* [Installation](#installation)
  * [Tools dependancies](#tools-dependancies)
    * [Use PEAR as read pairs merging software in preprocess](#use-pear-as-read-pairs-merging-software-in-preprocess)
  * [FROGS and dependancies installation](#frogs-and-dependancies-installation)
    * [From bioconda](#from-bioconda)
    * [From source](#from-source)
  * [Check intallation](#check-intallation)
* [Memory and parallelisation advices](#memory-and-parallelisation-advices)
* [Download databanks](#download-databanks)
* [Troubleshooting](#troubleshooting)
  * [Abnormal increase memory consumption with CPU number](#abnormal-increase-memory-consumption-with-cpu-number)
  * [Abnormal threads consumption in RDPClassifier](#abnormal-threads-consumption-in-rdpclassifier)
* [License](#license)
* [Copyright](#copyright)
* [Citation](#citation)
* [Contact](#contact)

# Convenient input data

Legend for the next schemas:
```
.: Complete nucleic sequence
!: Region of interest
*: PCR primers
```
* Paired-end classical protocol:
In the paired-end protocol R1 and R2 must share a nucleic region. 
For example the amplicons on 16S V3-V4 regions can have a length between 350 and 500nt, with 2*300pb sequencing the overlap is between 250nt and 100nt. 
```
        From:                                    To:
         rDNA .........!!!!!!................    ......!!!!!!!!!!!!!!!!!!!.....
         Ampl      ****!!!!!!****                  ****!!!!!!!!!!!!!!!!!!!****
           R1      --------------                  --------------
           R2      --------------                               --------------
```



In any case, the maximum overlap between R1 and R2 can be the complete overlap.

The minimum authorized overlap between R1 and R2 is 10nt. With less, the overlap can be incorrect, it will be rejected or considered as non overlap reads.

* Single-end classical protocol:

```
        rDNA .........!!!!!!................
        Ampl      ****!!!!!!****
        Read      --------------

```

* Custom protocol
```
        rDNA .....!!!!!!!!!!!!!!............
        Ampl      ****!!!!!!****
        Read      --------------       
```

The amplicons can have a high length variability such as ITS.  The R1 and R2 can have different length.

# Installation

This FROGS repository is for command line user. If you want to install FROGS on Galaxy, please refer to [FROGS-wrappers](https://github.com/geraldinepascal/FROGS-wrappers).

## Tools dependancies

FROGS is written in Python 3.7, uses home-made scripts written in PERL5 and R 3.6 and uses external Python library, numpy and Scipy.

FROGS relies on different specific tools for each of the analysis step.

| FROGS Tools                   |                          Dependancy                          |       version tested | last version |
| ----------------------------- | :----------------------------------------------------------: | -------------------: | ------------ |
| Preprocess and Remove_chimera |        [vsearch](https://github.com/torognes/vsearch)        | from 2.9.1 to 2.13.1 | 2.14.2       |
| Preprocess                    | [flash](https://sourceforge.net/projects/flashpage/files/) (optional) |               1.2.11 | last         |
| Preprocess                    |       [cutadapt](https://github.com/marcelm/cutadapt)        |                 1.18 | 2.9          |
| Clustering                    |          [swarm](https://github.com/torognes/swarm)          |                2.2.2 | 3.0.0        |
| ITSx                          |        [ITSx](http://microbiology.se/software/itsx/)         |      1.0.11 and 1.1b | last         |
| Affiliation_OTU               | [NCBI BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) |                2.7.1 | 2.9.0        |
| Affiliation_OTU               |    [RDP Classifier](https://github.com/rdpstaff/RDPTools)    |                2.0.3 | last         |
| Affiliation_OTU               | [EMBOSS needleall](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html) |                6.6.0 | last         |
| Tree                          |       [MAFFT](https://mafft.cbrc.jp/alignment/server/)       |                7.407 | 7.464        |
| Tree                          |     [Fasttree](http://www.microbesonline.org/fasttree/)      |               2.1.10 | last         |
| Tree / FROGSSTAT              | [plotly](https://plotly.com/r/), [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html), [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), [phyloseq](https://joey711.github.io/phyloseq/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [calibrate](https://cran.r-project.org/web/packages/calibrate/index.html), [formattable](https://cran.r-project.org/web/packages/formattable/), [DT](https://cran.r-project.org/web/packages/DT/index.html) |  depend on R version |              |

### **Use PEAR as read pairs merging software in preprocess**

[PEAR](https://cme.h-its.org/exelixis/web/software/pear/) is one of the most effective software for read pairs merging, but as its licence is not free for private use, we can not distribute it in FROGS.
If you work in an academic lab on a private Galaxy server, or if you have paid your licence you can use PEAR in FROGS preprocess.
For that you need to:

- have PEAR in your PATH or in the FROGS libexec directory. We have tested PEAR 0.9.10 version.
- use `--merge-software pear` option in the preprocess.py command line

## FROGS and dependancies installation

### From bioconda

FROGS is now available on bioconda (https://anaconda.org/bioconda/frogs).

  * to create a specific environment for a specific FROGS version

```
conda create --name __frogs@3.2.0 frogs=3.2.0
# to use FROGS, first you need to activate your environment
source activate __frogs@3.2.0
```

### From source

see [INSTALL_from_source.md](INSTALL_from_source.md)

## Check intallation

To check your installation you can type:
```
cd <FROGS_PATH>/test

sh test.sh ~/FROGS <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
# Note: JAVA_MEM must be at least 4 (= 4Gb of RAM).
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of Chaillou et al, ISME 2014, doi:10.1038/ismej.2014.202

This test executes the FROGS tools in command line mode.
Example:
```
[user@computer:/home/frogs/FROGS/test/]$ sh test.sh ~/FROGS 1 2 res
Step preprocess : Flash lundi 11 mai 2020, 17:31:35 (UTC+0200)
Step preprocess : Vsearch lundi 11 mai 2020, 17:34:48 (UTC+0200)
Step clustering lundi 11 mai 2020, 17:39:00 (UTC+0200)
Step remove_chimera lundi 11 mai 2020, 17:44:21 (UTC+0200)
Step otu filters lundi 11 mai 2020, 17:47:17 (UTC+0200)
Step ITSx lundi 11 mai 2020, 17:47:29 (UTC+0200)
Step affiliation_OTU lundi 11 mai 2020, 17:47:30 (UTC+0200)
Step affiliation_filter: masking mode lundi 11 mai 2020, 17:48:13 (UTC+0200)
Step affiliation_filter: deleted mode lundi 11 mai 2020, 17:48:13 (UTC+0200)
Step affiliation_postprocess lundi 11 mai 2020, 17:48:13 (UTC+0200)
Step clusters_stat lundi 11 mai 2020, 17:48:14 (UTC+0200)
Step affiliations_stat lundi 11 mai 2020, 17:48:16 (UTC+0200)
Step biom_to_tsv lundi 11 mai 2020, 17:48:22 (UTC+0200)
Step biom_to_stdBiom lundi 11 mai 2020, 17:48:22 (UTC+0200)
Step tsv_to_biom lundi 11 mai 2020, 17:48:22 (UTC+0200)
Step tree : mafft lundi 11 mai 2020, 17:48:23 (UTC+0200)
Step r_import_data lundi 11 mai 2020, 17:49:03 (UTC+0200)
Step r_composition lundi 11 mai 2020, 17:49:50 (UTC+0200)
Step r_alpha_diversity lundi 11 mai 2020, 17:50:37 (UTC+0200)
Step r_beta_diversity lundi 11 mai 2020, 17:51:05 (UTC+0200)
Step r_structure lundi 11 mai 2020, 17:51:21 (UTC+0200)
Step r_clustering lundi 11 mai 2020, 17:52:06 (UTC+0200)
Step r_manova lundi 11 mai 2020, 17:52:21 (UTC+0200)
Step deseq2_preprocess lundi 11 mai 2020, 17:52:35 (UTC+0200)
Step deseq2_visualization lundi 11 mai 2020, 17:53:06 (UTC+0200)
Completed with success
```



# Memory and parallelisation advices

If you have more than one CPU, it is recommended to increase the number of CPUs used by tools.
All the CPUs must be on the same computer/node.

|         Tool          | RAM per CPU | Minimal RAM | Configuration example |
| :-------------------: | :---------: | :---------: | :-------------------: |
|      Preprocess       |     8Gb     |      -      |   12 CPUs and 96 GB   |
|      Clustering       |      -      |    10 Gb    |   16 CPUs and 60 GB   |
| ITSx / Remove_Chimera |     3Gb     |     5Gb     |   12 CPUs and 36 GB   |
|    Affiliation_OTU    |      -      |    20 Gb    |  30 CPUs and 300 GB   |



# Download databanks

Reference database are needed to filter contaminants, assign taxonomy to each OTU or filter ambiguities for hyper variable amplicon length.

We propose some databanks, that you simply need to download and extract.

Please take time to read individual README.txt and LICENCE.txt files.

* Assignation databank

  these databanks are formatted for NCBI Blast+ and RDP Classifier

  [available databases](http://genoweb.toulouse.inrae.fr/frogs_databanks/assignation/readme.txt) : http://genoweb.toulouse.inra.fr/frogs_databanks/assignation

* Contaminant databank

  these banks are formatted for NCBI Blast+

  http://genoweb.toulouse.inrae.fr/frogs_databanks/contaminants

* Hyper variable in length amplicon databank

  This is simply fasta file.

  http://genoweb.toulouse.inrae.fr/frogs_databanks/HVL




# Troubleshooting
## Abnormal increase memory consumption with CPU number
With some old versions of glibc the virtual memory used by CPU is multiplicative.

| Nb CPUs | expected RAM consumtion | observed RAM consumption |
| :-----: | :---------------------: | :----------------------: |
|    1    |          1 Gb           |           1Gb            |
|    2    |          2 Gb           |          2*2 Gb          |
|    3    |          3 Gb           |          3*3 Gb          |
|    4    |          4 Gb           |          4*4 Gb          |


The parameters memory and CPU provided in examples take into account this problem.

## Abnormal threads consumption in RDPClassifier
With large databases like silva NR the RDPClassifier opens automatically a large number of threads. These threads use all the available CPU ressources. This is not an acceptable behaviour in multi-user context.
To prevent this behaviour the tool 'affiliation_OTU' uses taskset to force RDPClassifier to run only on the specified number of CPUs. The number of threads is not changed but the CPU consumption is controled.


# License
GNU GPL v3


# Copyright
2020 INRAE


# Citation
Please cite the **FROGS article**: *Escudie F., et al. Bioinformatics, 2018. FROGS: Find, Rapidly, OTUs with Galaxy Solution.*


# Contact
frogs@inrae.fr
