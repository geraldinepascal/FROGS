<p align="center">
 <a href="http://frogs.toulouse.inra.fr/">
  <img src="FROGS_logo.png" align="center" width="20%" style="display: block; margin: auto;"/>
 </a>
</p>

Visit our web site : http://frogs.toulouse.inrae.fr/

[![Release](https://img.shields.io/badge/release-3.2.2-blue.svg)![Date](https://img.shields.io/badge/date-March%202021-red.svg)](https://github.com/geraldinepascal/FROGS-wrappers/releases) [<img src="https://www.podcastscience.fm/wp-content/uploads/2017/12/deezer.png" width="5%" style="display: block; margin: auto;"/>](https://www.deezer.com/fr/playlist/5233843102?utm_source=deezer&utm_content=playlist-5233843102&utm_term=18632989_1545296531&utm_medium=web)



# Description

FROGS is a CLI workflow designed to produce an OTU count matrix from high depth sequencing amplicon data.

FROGS-wrappers allow to add FROGS on a Galaxy instance. (see https://github.com/geraldinepascal/FROGS-wrappers)

This workflow is focused on:

- User-friendliness with lots of rich graphic outputs and the integration in Galaxy thanks to FROGS-wrappers.
- Accuracy with a clustering without global similarity threshold, the management of separated PCRs in the chimera removal step, and the management of multi-affiliations.
- Dealing of non overlapping pair of sequences from long amplicon like ITS, or RPB2.
- Speed with fast algorithms parallelisation and easy to use.
- Scalability with algorithms designed to support the data growth.



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
In the paired-end protocol R1 and R2 may share a nucleic region. 
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

FROGS is written in Python 3.7 (with external numpy and Scipy libraries) , uses also home-made scripts written in PERL5 and R 3.6.

FROGS relies on different specific tools for each of the analysis steps.

| FROGS Tools |Dependancy  | version tested |
| ----------- | :--------: | -------------: |
| Preprocess and Remove_chimera |        [vsearch](https://github.com/torognes/vsearch)        | 2.15.1 |
| Preprocess                    | [flash](https://sourceforge.net/projects/flashpage/files/) (optional) |               1.2.11 |
| Preprocess                    |       [cutadapt](https://github.com/marcelm/cutadapt)        |            3.1 |
| Clustering                    |          [swarm](https://github.com/torognes/swarm)          |            3.0.0 |
| ITSx                          |        [ITSx](http://microbiology.se/software/itsx/)         |  1.1.2 |
| Affiliation_OTU               | [NCBI BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) |          2.10.1 |
| Affiliation_OTU               |    [RDP Classifier](https://github.com/rdpstaff/RDPTools)    |                2.0.3 |
| Affiliation_OTU               | [EMBOSS needleall](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html) |                6.6.0 |
| Tree                          |       [MAFFT](https://mafft.cbrc.jp/alignment/server/)       |                7.475 |
| Tree                          |     [Fasttree](http://www.microbesonline.org/fasttree/)      |               2.1.10 |
| Tree / FROGSSTAT              | [plotly](https://plotly.com/r/), [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html), [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), [phyloseq](https://joey711.github.io/phyloseq/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [calibrate](https://cran.r-project.org/web/packages/calibrate/index.html), [formattable](https://cran.r-project.org/web/packages/formattable/), [DT](https://cran.r-project.org/web/packages/DT/index.html) |              R 3.6.3 |
| FROGSSTAT | [pandoc](https://pandoc.org/) | 2.11.3|

### **Use PEAR as read pairs merging software in preprocess**

[PEAR](https://cme.h-its.org/exelixis/web/software/pear/) is one of the most effective software for read pairs merging, but as its license is not free for private use, we can not distribute it in FROGS.
If you work in an academic lab on a private Galaxy server, or if you have paid your license you can use PEAR in FROGS preprocess.
For that you need to:

- have PEAR in your PATH or in the FROGS libexec directory. We have tested PEAR 0.9.10 version (last version 0.9.11).
- use `--merge-software pear` option in the preprocess.py command line

## FROGS and dependancies installation

### From conda

FROGS is now available on bioconda (https://anaconda.org/bioconda/frogs).

  * to create a specific environment for a specific FROGS version

```
conda env create --name frogs@3.2.2 --file frogs-conda-requirements.yaml
# to use FROGS, first you need to activate your environment
conda activate frogs@3.2.2
```

### From source

see [INSTALL_from_source.md](INSTALL_from_source.md)

## Check intallation

To check your installation you can type:
```
cd <FROGS_PATH>/test
# when using conda FROGS_PATH=<conda_env_dir>/__frogs@3.2.2/share/FROGS_3.2.2

sh test.sh <FROGS_PATH> <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of [Chaillou et al, ISME 2014](https://doi.org/10.1038/ismej.2014.202)

This test executes the FROGS tools in command line mode.
Example:

```
[user@computer:/home/frogs/FROGS/test/]$ sh test.sh ../ 1 2 res
Step preprocess : Flash mardi 10 novembre 2020, 10:56:56 (UTC+0100)
Step preprocess : Vsearch mardi 10 novembre 2020, 10:59:57 (UTC+0100)
Step clustering mardi 10 novembre 2020, 11:02:51 (UTC+0100)
Step remove_chimera mardi 10 novembre 2020, 11:08:31 (UTC+0100)
Step otu filters mardi 10 novembre 2020, 11:13:43 (UTC+0100)
Step ITSx mardi 10 novembre 2020, 11:14:00 (UTC+0100)
Step affiliation_OTU mardi 10 novembre 2020, 11:14:01 (UTC+0100)
Step affiliation_filter: masking mode mardi 10 novembre 2020, 11:14:53 (UTC+0100)
Step affiliation_filter: deleted mode mardi 10 novembre 2020, 11:14:54 (UTC+0100)
Step affiliation_postprocess mardi 10 novembre 2020, 11:14:54 (UTC+0100)
Step normalisation mardi 10 novembre 2020, 11:14:55 (UTC+0100)
Step clusters_stat mardi 10 novembre 2020, 11:14:55 (UTC+0100)
Step affiliations_stat mardi 10 novembre 2020, 11:14:58 (UTC+0100)
Step biom_to_tsv mardi 10 novembre 2020, 11:15:05 (UTC+0100)
Step biom_to_stdBiom mardi 10 novembre 2020, 11:15:06 (UTC+0100)
Step tsv_to_biom mardi 10 novembre 2020, 11:15:06 (UTC+0100)
Step tree mardi 10 novembre 2020, 11:15:06 (UTC+0100)
Step phyloseq_import_data mardi 10 novembre 2020, 11:16:36 (UTC+0100)
Step phyloseq_composition mardi 10 novembre 2020, 11:18:00 (UTC+0100)
Step phyloseq_alpha_diversity mardi 10 novembre 2020, 11:19:31 (UTC+0100)
Step phyloseq_beta_diversity mardi 10 novembre 2020, 11:20:19 (UTC+0100)
Step phyloseq_structure mardi 10 novembre 2020, 11:20:45 (UTC+0100)
Step phyloseq_clustering mardi 10 novembre 2020, 11:21:59 (UTC+0100)
Step phyloseq_manova mardi 10 novembre 2020, 11:22:20 (UTC+0100)
Step deseq2_preprocess mardi 10 novembre 2020, 11:22:42 (UTC+0100)
Step deseq2_visualisation mardi 10 novembre 2020, 11:23:29 (UTC+0100)
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


# License
GNU GPL v3


# Copyright
2020 INRAE


# Citation
Please cite the [**FROGS article**: *Escudie F., et al. Bioinformatics, 2018. FROGS: Find, Rapidly, OTUs with Galaxy Solution.*](https://doi.org/10.1093/bioinformatics/btx791)


# Contact
frogs-support@inrae.fr
