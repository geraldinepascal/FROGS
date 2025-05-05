<p align="center">
 <a href="http://frogs.toulouse.inrae.fr/">
  <img src="FROGS_logo.png" align="center" width="20%" style="display: block; margin: auto;"/>
 </a>
</p>

Visit our web site : http://frogs.toulouse.inrae.fr/

[![Release](https://img.shields.io/badge/release-5.0.2-blue.svg)![Date](https://img.shields.io/badge/date-May%202025-red.svg)](https://github.com/geraldinepascal/FROGS-wrappers/releases) [<img src="https://cdn.icon-icons.com/icons2/70/PNG/512/deezer_14086.png" width="2%" style="display: block; margin: auto;"/>](https://www.deezer.com/fr/playlist/5233843102?utm_source=deezer&utm_content=playlist-5233843102&utm_term=18632989_1545296531&utm_medium=web)



# Description

FROGS is a CLI workflow designed to produce and analysis an ASV count matrix from high depth sequencing amplicon data.

FROGS-wrappers allow to add FROGS on a Galaxy instance. (see https://github.com/geraldinepascal/FROGS-wrappers)

This workflow is focused on:

- User-friendliness with lots of rich graphic outputs and the integration in Galaxy thanks to FROGS-wrappers.
- Dealing with **Illumina**, **IonTorrent**, **454**, **PacBio** and **Oxford Nanopore** sequencing technology reads.
- Dealing of non overlapping pair of sequences from long amplicon like ITS, or RPB2.
- Accuracy with a clustering without global similarity threshold or a denoising process, the management of separated PCRs in the chimera removal step, and the management of multi-affiliations.
- Access to general statistics on microbial communities and differential abundance analysis.
- Prediction of functionnal annotation and metabolic pathways.
- Speed with fast algorithms parallelisation and easy to use.
- Scalability with algorithms designed to support the data growth.

# Table of content

* [Convenient input data](#convenient-input-data)
* [Installation](#installation)
  * [Tools dependencies](#tools-dependencies)
    * [Use PEAR as read pairs merging software in preprocess](#use-pear-as-read-pairs-merging-software-in-preprocess)
  * [FROGS and dependencies installation](#frogs-and-dependencies-installation)
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

## Tools dependencies

FROGS is written in Python 3 (with external numpy and Scipy libraries) , uses also home-made scripts written in PERL5 and R 4.

FROGS relies on different specific tools for each of the analysis steps.

| FROGS Tools |Dependancy  | version tested |
| ----------- | :--------: | -------------: |
| Denoising and Remove_chimera |        [vsearch](https://github.com/torognes/vsearch)        | 2.17.0 |
| Denoising                    | [flash](https://sourceforge.net/projects/flashpage/files/) (optional) |               1.2.11 |
| Denoising                    |       [cutadapt](https://github.com/marcelm/cutadapt) (need to be >=2.8)       |            2.10 |
| Denoising                    |          [swarm](https://github.com/torognes/swarm) (need to be >=2.1)          |            3.1.4 |
| Denoising                    |          [DADA2](https://benjjneb.github.io/dada2/index.html)     |            1.22.0 |
| ITSx                          |        [ITSx](http://microbiology.se/software/itsx/)         |  1.1.2 |
| Taxonomic_affiliation               | [NCBI BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) |          2.16 |
| Taxonomic_affiliation               | [EMBOSS needleall](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html) |                6.6.0 |
| Tree                          |       [MAFFT](https://mafft.cbrc.jp/alignment/server/)       |                7.525 |
| Tree                          |     [Fasttree](http://www.microbesonline.org/fasttree/)      |               2.1.9 |
| Tree / FROGSSTAT              | [plotly](https://plotly.com/r/), [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html), [rmarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), [phyloseq](https://joey711.github.io/phyloseq/), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [calibrate](https://cran.r-project.org/web/packages/calibrate/index.html), [formattable](https://cran.r-project.org/web/packages/formattable/), [DT](https://cran.r-project.org/web/packages/DT/index.html) |              R 4.1.2 |
| FROGSSTAT | [pandoc](https://pandoc.org/) | 2.11.3|
| FROGSFUNC | [PICRUSt2](https://github.com/picrust/picrust2/wiki) | 2.5.1|
| FROGSFUNC | [ete3](http://etetoolkit.org/) | 3.1.1|

### **Use PEAR as read pairs merging software in preprocess**

[PEAR](https://cme.h-its.org/exelixis/web/software/pear/) is one of the most effective software for read pairs merging, but as its license is not free for private use, we can not distribute it in FROGS.
If you work in an academic lab on a private Galaxy server, or if you have paid your license you can use PEAR in FROGS preprocess.
For that you need to:

- have PEAR in your PATH or in the FROGS libexec directory. We have tested PEAR 0.9.11 version.
- use `--merge-software pear` option in the preprocess.py command line

## FROGS and dependencies installation

### From conda

FROGS is now available on bioconda (https://anaconda.org/bioconda/frogs).

  * to create a specific environment for a specific FROGS version

```
conda env create --name frogs@5.0.2 --file frogs-conda-requirements.yaml
# to use FROGS, first you need to activate your environment
conda activate frogs@5.0.2
```

**WARNING** : As PICRUSt2 currently relies on a different R version, in order to use the FROGSFUNC tools, it is necessary to create a dedicated conda environment as follows:

```
conda env create --name frogsfunc@5.0.2 --file frogsfunc-conda-requirements.yaml
# and then activate the environment
conda activate frogsfunc@5.0.2
```

After that, you just have to switch from one environment to another (with `conda activate frogs@5.0.2` or `conda activate frogsfunc@5.0.2` depending on whether you want to use FROGSFUNC or all the other tools.

## Check intallation

To check your installation you can type:

```
cd <conda_env_dir>/frogs@5.0.2/share/FROGS-5.0.2/test

conda activate frogs@5.0.2

sh test_frogs.sh <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
```
"Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
"Statistical" tools are performed on an extract of the published results of [Chaillou et al, ISME 2014](https://doi.org/10.1038/ismej.2014.202)

This test executes the FROGS tools in command line mode.
Example:

```
[user@computer:/home/frogs/FROGS/test/]$ sh test_frogs.sh 1 2 res
Step demultiplex Wed Apr 30 11:40:08 AM CEST 2025
Step denoising 16S vsearch Wed Apr 30 11:40:09 AM CEST 2025:
Step denoising 16S vsearch swarm denoising and distance 3 Wed Apr 30 11:40:25 AM CEST 2025:
Step denoising 16S pear Wed Apr 30 11:40:41 AM CEST 2025:
Step denoising: dada2 keep-unmerged Wed Apr 30 11:41:49 AM CEST 2025
Step denoising: preprocess only Wed Apr 30 11:42:12 AM CEST 2025
Step remove_chimera Wed Apr 30 11:42:16 AM CEST 2025
Step cluster_filters Wed Apr 30 11:42:19 AM CEST 2025
Step itsx Wed Apr 30 11:42:26 AM CEST 2025
Step taxonomic_affiliation Wed Apr 30 11:42:31 AM CEST 2025
Step affiliation_filters: masking mode Wed Apr 30 11:42:43 AM CEST 2025
Step affiliation_filters: deleted mode Wed Apr 30 11:42:47 AM CEST 2025
Step affiliation_postprocess Wed Apr 30 11:42:50 AM CEST 2025
Step normalisation Wed Apr 30 11:42:52 AM CEST 2025
Step cluster_stats Wed Apr 30 11:42:57 AM CEST 2025
Step affiliation_stats Wed Apr 30 11:42:58 AM CEST 2025
Step biom_to_tsv Wed Apr 30 11:42:59 AM CEST 2025
Step biom_to_stdBiom Wed Apr 30 11:43:00 AM CEST 2025
Step tsv_to_biom Wed Apr 30 11:43:00 AM CEST 2025
Step tree Wed Apr 30 11:43:01 AM CEST 2025
Step phyloseq_import_data Wed Apr 30 11:43:13 AM CEST 2025
Step phyloseq_composition Wed Apr 30 11:44:20 AM CEST 2025
Step phyloseq_alpha_diversity Wed Apr 30 11:44:44 AM CEST 2025
Step phyloseq_beta_diversity Wed Apr 30 11:45:08 AM CEST 2025
Step phyloseq_structure Wed Apr 30 11:45:24 AM CEST 2025
Step phyloseq_clustering Wed Apr 30 11:45:39 AM CEST 2025
Step phyloseq_manova Wed Apr 30 11:45:53 AM CEST 2025
Step deseq2_preprocess Wed Apr 30 11:46:08 AM CEST 2025
DESeq2 asv abundances
DESeq2 function abundances
Step deseq2_visualisation Wed Apr 30 11:47:30 AM CEST 2025
DESeq2 ASV abundances
DESeq2 function abundances
Completed with success
```

Finally, to check the FROGSFUNC tools installation you can type:

```
cd <conda_env_dir>/frogsfunc@5.0.2/share/FROGS-5.0.2/test

conda activate frogsfunc@5.0.2

sh test_frogsfunc.sh <OUT_FOLDER>
```

This test executes the FROGSFUNC tools in command line mode.
Example:

```
[user@computer:/home/frogs/FROGS/test/]$ sh test_frogsfunc.sh res
Step frogsfunc_placeseqs Wed Apr 30 11:50:43 AM CEST 2025
Step frogsfunc_functions Wed Apr 30 11:52:04 AM CEST 2025
Step frogsfunc_pathways Wed Apr 30 11:54:47 AM CEST 2025
Completed with success
```


# Memory and parallelisation advices

If you have more than one CPU, it is recommended to increase the number of CPUs used by tools.
All the CPUs must be on the same computer/node.

|         Tool          | RAM per CPU | Minimal RAM | Configuration example |
| :-------------------: | :---------: | :---------: | :-------------------: |
|      Denoising       |     8Gb     |      10 Gb      |   12 CPUs and 96 GB   |
| ITSx / Remove_Chimera |     3Gb     |     5Gb     |   12 CPUs and 36 GB   |
|    Taxonomic_affiliation    |      -      |    20 Gb    |  30 CPUs and 300 GB   |



# Download databanks

Reference databanks are needed to filter contaminants, assign taxonomy to each ASV or filter ambiguities for hyper variable amplicon length.

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

In addition, several default databases are used in FROGSFUNC steps by PICRUSt2.


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
2025 INRAE


# Citation

Depending on which type of amplicon you are working on (mergeable or unmergeable), please cite  one of the two FROGS publications:
* [*Escudie F., et al. Bioinformatics, 2018. FROGS: Find, Rapidly, OTUs with Galaxy Solution.*](https://doi.org/10.1093/bioinformatics/btx791)
* [*Bernard M., et al. Briefings in Bioinformatics, 2021. FROGS: a powerful tool to analyse the diversity of fungi with special management of internal transcribed spacers.*](https://doi.org/10.1093/bib/bbab318)

# Contact
frogs-support@inrae.fr
