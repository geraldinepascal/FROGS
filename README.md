[<img src="FROGS_logo.png" width="20%" style="display: block; margin: auto;"/>](http://frogs.toulouse.inra.fr/)

Visit our web site : http://frogs.toulouse.inra.fr/



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

FROGS is written in Python 2.7, uses home-made scripts written in PERL5 and R3.6 and uses external Python library, numpy and Scipy.

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

  [available databases](http://genoweb.toulouse.inra.fr/frogs_databanks/assignation/readme.txt) : http://genoweb.toulouse.inra.fr/frogs_databanks/assignation

* Contaminant databank

  these banks are formatted for NCBI Blast+

  http://genoweb.toulouse.inra.fr/frogs_databanks/contaminants

* Hyper variable in length amplicon databank

  This is simply fasta file.

  http://genoweb.toulouse.inra.fr/frogs_databanks/HVL




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
2018 INRA


# Citation
Please cite the **FROGS article**: *Escudie F., et al. Bioinformatics, 2018. FROGS: Find, Rapidly, OTUs with Galaxy Solution.*


<<<<<<< HEAD
    Note: The amplicons can have a length variability.
          The R1 and R2 can have different length.


## Installation
### 1. Download code
    Released versions
        Available at https://github.com/geraldinepascal/FROGS/releases
        After the download, follow the README instructions.

    Otherwise, you can clone the repository from GitHub:
        git clone https://github.com/geraldinepascal/FROGS.git

### 2. FROGS folder
    Default structure:
        <FROGS_PATH>/
            app/                 # FROGS applications (this folder must be 
                ...              # accessible for command line and/or galaxy)
                preprocess.py        # Link to tools/preprocess/preprocess.py
                ...
            lib/                 # FROGS python librairies
                ...
                frogsBiom.py
                ...
            libexec/             # FROGS softwares (you can also add the 
                ...              # dependencies in this folder)
                biomTools.py
                ...
            tools/               # FROGS applications with one sub-folder by
                ...              # application
                preprocess/
                    preprocess.py
                    preprocess.xml
                ...

    If you want to change this architecture 'libexec' must be accessible in the
    PATH and 'lib' must be accessible in the PYTHONPATH.

### 3. Install dependencies

### 3.1 lib exec and additionnal packages

    Dependencies must be accessible in the PATH or added in <FROGS_PATH>/libexec.
    
    python interpreter
        Version: 2.7
        Tools: all

    python SciPy
        Tools: clusters_stat

    perl interpreter
        Version: 5
        Tools: demultiplex
    
    vsearch
        Version: 2.6.2
        Named as: vsearch
        Tools: preprocess and remove_chimera
        Download: https://github.com/torognes/vsearch
        Warning : zlib and bzlib need to be installed before compilling vsearch to deal with fastq.gz or fastq.bz2 files.

    flash (optional, but recommended for user that used FROGS 2.0)
        Version: 1.2.11
        Named as: flash
        Tools: preprocess
        Download: https://sourceforge.net/projects/flashpage/files/

    pear (optional)
        Version: 0.9.10
        Named as: pear
        Tools: preprocess
        Download: https://sco.h-its.org/exelixis/web/software/pear/

    cutadapt
        Version: 1.8.3
            Note : With the cutadapt version 1.12, the memory usage increases drastically. 
                   We advise our user to install, at most, the cutadapt version  1.11       
        Named as: cutadapt
        Tools: preprocess
        Download: https://github.com/marcelm/cutadapt
                  OR
                  https://pypi.python.org/pypi/cutadapt
        
    swarm
        Version: 2.1.1
        Named as: swarm
        Tools: clustering
        Download: https://github.com/torognes/swarm

    ITSx
        Version : 1.0.11
        Named : ITSx
        Tools : itsx
        Download : http://microbiology.se/software/itsx/
        Remark : ITSx_db folder need to be in the PATH or in <FROGS_PATH>/libexec
                 it depends on HMMER 3 or later (only for hmmpress and hmmscan need to be linked in <FROGS_PATH>/libexec or available in the PATH)
                 if ITSx test command line failed it's may be due to a difference in HMMER version used to prepare HMM models: 
                        cd <ITSx_DIR> ; rm ITSx_db/HMMs/*.h3*
                        for hmm in `ls  ITSx_db/HMMs/*.hmm `; do hmmpress $hmm; done

    NCBI Blast+ blastn
        Version: 2.2.30+
        Named as: blastn
        Tools: affiliation_OTU and filters
        Download: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

    RDPClassifier
        Version: 2.0.2
        Named as: classifier.jar
        Tools: affiliation_OTU
        Download: https://github.com/rdpstaff/RDPTools

    taskset
        Version: util-linux-ng 2.17.2
        Named as: taskset
        Tools: affiliation_OTU
        Install: sudo apt-get install util-linux
                 OR
                 sudo yum install util-linux

    Needlall
        Version: EMBOSS:6.6.0.0
        Named as: needleall
        Tools : affiliation_OTU
        Download : ftp://emboss.open-bio.org/pub/EMBOSS/

    Pynast
        Version: 1.2.2
        Named as: pynast
        Tools: tree
        Download: https://pypi.python.org/pypi/pynast

    Mafft
        Version: v7.310
        Named as: mafft
        Tools: tree
        Download: http://mafft.cbrc.jp/alignment/software/

    Fasttree
        Version: 2.1.10
        Named as: FastTree
        Tools: tree
        Download: http://www.microbesonline.org/fasttree/#Install

    RScript
        Version : 3.4.0
        Named as : RScript
        Tools : all FROGSSTAT_Phyloseq
        Download : https://cran.r-project.org/

    Phangorn R package
        Version : 2.2.0
        Tools : FROGS_Tree
        Installation in R session : # https://cran.r-project.org/web/packages/phangorn/index.html
                                    install.packages("phangorn")
        Test in R session : library(phangorn)

    Rmarkdown R package
        Version : 1.5
        Tools : all FROGSSTAT_Phyloseq
        Install in R session : # https://cran.r-project.org/web/packages/rmarkdown/index.html
                               install.packages("rmarkdown")

    Pandoc
        Version : 1.17.2
        Named as : pandoc
        Tools : all FROGSSTAT_Phyloseq
        Download/Installation : # http://pandoc.org/installing.html#linux or simply soft-link pandoc binary from RStudio path (if you have Rstudio installed)

    Phyloseq R package
        Version : 1.20.0
        Tools : all FROGSSTAT_Phyloseq
        Installation in R session : # https://joey711.github.io/phyloseq/install.html
                                    source("https://bioconductor.org/biocLite.R") ; biocLite("phyloseq")
        Test in R session : library(phyloseq)

    Plotly R package
        Version : 4.7.0
        Tools : FROGSSTAT_Phyloseq_composition and FROGSSTAT_Phyloseq_structure
        Installation in R session : # https://plot.ly/r/getting-started/
                                    install.packages("plotly")
        Test in R session : library(plotly)

    GridExtra R package
        Version : 2.2.1
        Tools : FROGSSTAT_Phyloseq_Beta_Diversity, FROGSSTAT_Phyloseq_Sample_Clustering, FROGSSTAT_Phyloseq_composition and FROGSSTAT_Phyloseq_structure
        Installation in R session : # https://cran.r-project.org/web/packages/gridExtra/index.html
                                    install.packages("gridExtra")
        Test in R session : library(gridExtra)

### 3.2 R lib 

    Dependencies must be accessible in <FROGS_PATH>/lib/external-lib.

    Phyloseq-extended
        Version : v0.99
        Tools : all FROGSSTAT tools
        Installation : # https://github.com/mahendra-mariadassou/phyloseq-extended/releases
                       untar archive and copy or link content of folder "phyloseq-extended/" in <FROGS_PATH>/lib/external-lib


### 4. Check intallation
    To check your installation you can type:
        cd <FROGS_PATH>/test
        bash test.sh <FROGS_PATH> <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
    
    "Bioinformatic" tools are performed on a small simulated dataset of one sample replicated three times.
    "Statistical" tools are performed on an extract of the published results of Chaillou et al, ISME 2014, doi:10.1038/ismej.2014.202

    This test executes the FROGS tools in command line mode.
    Note:
        JAVA_MEM must be at least 4 (= 4Gb of RAM).
    Example:
        [user@computer:/home/user]$cd /home/user/frogs_git/test
        [user@computer:/home/user/frogs_git/test]$bash test.sh /home/user/frogs_git/ 2 4 /tmp/results
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
        Step tree : pynast mercredi 10 octobre 2018, 14:23:12 (UTC+0200)
        Step tree : mafft mercredi 10 octobre 2018, 14:23:26 (UTC+0200)
        Step r_import_data mercredi 10 octobre 2018, 14:25:25 (UTC+0200)
        Step r_composition mercredi 10 octobre 2018, 14:25:39 (UTC+0200)
        Step r_alpha_diversity mercredi 10 octobre 2018, 14:25:53 (UTC+0200)
        Step r_beta_diversity mercredi 10 octobre 2018, 14:26:19 (UTC+0200)
        Step r_structure mercredi 10 octobre 2018, 14:26:31 (UTC+0200)
        Step r_clustering mercredi 10 octobre 2018, 14:26:47 (UTC+0200)
        Step r_manova mercredi 10 octobre 2018, 14:26:57 (UTC+0200)
        Completed with success

### 5. New datatype
    
#### 5.1 Add the biom1 datatype in galaxy
    The biom1 datatype is availaible since version 16 of galaxy.

    For previous versions add the following text in galaxy datatypes_conf.xml:
        For galaxy releases 15:
            <registration converters_path="lib/galaxy/datatypes/converters">
                ...
                <datatype extension="biom1" type="galaxy.datatypes.text:Biom1" display_in_upload="True" subclass="True" mimetype="application/json" />
                ...
            <registration />
        For galaxy releases 14:
            <registration converters_path="lib/galaxy/datatypes/converters">
                ...
                <datatype extension="biom1" type="galaxy.datatypes.data:Text" display_in_upload="True" subclass="True" mimetype="application/json" />
                ...
            <registration />

#### 5.2 Add the tar datatype in galaxy
    Datatype tar is available since release 16.07 of galaxy.

    For previous 16.XX version add the following text in galaxy datatypes_conf.xml:
        <datatype extension="tar" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>

    For previous version you can use the FROGS_Upload_archive tool (this tool is no more working in version 16.XX and later).

### 6. Add tools in galaxy
    Add the tools in galaxy tool_conf.xml.
    Example:
        ...
        <section id="FROGS_wrappers" name="FROGS">
        <label text="OTUs reconstruction" id="FROGS_OTU" />
            <tool file="FROGS/app/upload_tar.xml" />
            <tool file="FROGS/app/demultiplex.xml" />
            <tool file="FROGS/app/preprocess.xml" />
            <tool file="FROGS/app/clustering.xml" />
            <tool file="FROGS/app/remove_chimera.xml" />  
            <tool file="FROGS/app/filters.xml" />
            <tool file="FROGS/app/itsx.xml" />
            <tool file="FROGS/app/affiliation_OTU.xml" />
            <tool file="FROGS/app/affiliation_postprocess.xml" />
            <tool file="FROGS/app/normalisation.xml" />
            <tool file="FROGS/app/clusters_stat.xml" />
            <tool file="FROGS/app/affiliations_stat.xml" />
            <tool file="FROGS/app/biom_to_stdBiom.xml" />
            <tool file="FROGS/app/biom_to_tsv.xml" />
            <tool file="FROGS/app/tsv_to_biom.xml" />
            <tool file="FROGS/app/tree.xml" />
        <label text="OTUs structure and composition analysis" id="FROGSSTAT_Phyloseq" />
            <tool file="FROGS/app/r_import_data.xml" />
            <tool file="FROGS/app/r_composition.xml" />
            <tool file="FROGS/app/r_alpha_diversity.xml" />
            <tool file="FROGS/app/r_beta_diversity.xml" />
            <tool file="FROGS/app/r_structure.xml" />
            <tool file="FROGS/app/r_clustering.xml" />
            <tool file="FROGS/app/r_manova.xml" />
        </section>
        ...
    Note: 
        <tool file="FROGS/app/upload_tar.xml" /> is no more working in galaxy version 16.XX and later. 
        Prefer to add tar datatype.
        
### 7. Set memory and parallelisation settings
    If you have more than one CPU, it is recommended to increase the number
    of CPUs used by tools.
    All the CPUs must be on the same computer/node.

    a] Specifications  
        Tool            RAM/CPU     Minimal RAM     Configuration example
        affiliation        -          20 Gb          30 CPUs and 300 GB
        chimera          3 Gb          5 Gb          12 CPUs and 36 GB
        clustering         -          10 Gb          16 CPUs and 60 GB
        preprocess       8 Gb            -           12 CPUs and 96 GB

    b] Change the tool launcher configuration
        In galaxy the job_conf.xml allows to change the scheduler 
        submission parameters.
        Example for SGE scheduler:
            <destinations>
                ...
                <destination id="FROGS_preprocess_job" runner="drmaa">
                    <param id="galaxy_external_runjob_script">scripts/drmaa_external_runner.py</param>
                    <param id="galaxy_external_killjob_script">scripts/drmaa_external_killer.py</param>
                    <param id="galaxy_external_chown_script">scripts/external_chown_script.py</param>
                    <param id="nativeSpecification">-clear -q galaxyq -l mem=5G -l h_vmem=13G -pe parallel_smp 12</param>
                </destination>
                <destination id="FROGS_clustering_job" runner="drmaa">
                    <param id="galaxy_external_runjob_script">scripts/drmaa_external_runner.py</param>
                    <param id="galaxy_external_killjob_script">scripts/drmaa_external_killer.py</param>
                    <param id="galaxy_external_chown_script">scripts/external_chown_script.py</param>
                    <param id="nativeSpecification">-clear -q galaxyq -l mem=3G -l h_vmem=10G -pe parallel_smp 16</param>
                </destination>
                <destination id="FROGS_remove_chimera_job" runner="drmaa">
                    <param id="galaxy_external_runjob_script">scripts/drmaa_external_runner.py</param>
                    <param id="galaxy_external_killjob_script">scripts/drmaa_external_killer.py</param>
                    <param id="galaxy_external_chown_script">scripts/external_chown_script.py</param>
                    <param id="nativeSpecification">-clear -q galaxyq -l mem=3G -l h_vmem=4G -pe parallel_smp 12</param>
                </destination>
                <destination id="FROGS_itsx_job" runner="drmaa">
                    <param id="galaxy_external_runjob_script">scripts/drmaa_external_runner.py</param>
                    <param id="galaxy_external_killjob_script">scripts/drmaa_external_killer.py</param>
                    <param id="galaxy_external_chown_script">scripts/external_chown_script.py</param>
                    <param id="nativeSpecification">-clear -q galaxyq -l mem=3G -l h_vmem=4G -pe parallel_smp 12</param>
                </destination>
                <destination id="FROGS_affiliation_OTU_job" runner="drmaa">
                    <param id="galaxy_external_runjob_script">scripts/drmaa_external_runner.py</param>
                    <param id="galaxy_external_killjob_script">scripts/drmaa_external_killer.py</param>
                    <param id="galaxy_external_chown_script">scripts/external_chown_script.py</param>
                    <param id="nativeSpecification">-clear -q galaxyq -l mem=7G -l h_vmem=10G -pe parallel_smp 30</param>
                </destination>
            </destinations>
            <tools>
                ...
                <tool id="FROGS_preprocess" destination="FROGS_preprocess_job"/>   
                <tool id="FROGS_clustering" destination="FROGS_clustering_job"/>     
                <tool id="FROGS_remove_chimera" destination="FROGS_remove_chimera_job"/> 
                <tool id="FROGS_itsx" destination="FROGS_itsx_job"/> 
                <tool id="FROGS_affiliation_OTU" destination="FROGS_affiliation_OTU_job"/>
            </tools>

### 8. Upload and configure the databanks
    a] Assignation databank
        - Upload databanks and indexes from http://genoweb.toulouse.inra.fr/frogs_databanks/assignation
        - Extract databanks.
        - To use these databank, you need to create a .loc file named
          'frogs_db.loc'. The path provided must be the '.fasta'.
          (see the frogs_db.loc example file)
          
    b] Contaminant databank
        - Upload databank and indexes from http://genoweb.toulouse.inra.fr/frogs_databanks/contaminants
        - Extract databank.
        - To use this databank, you need to create a .loc file named
          'phiX_db.loc'. The path provided must be the '.fasta'.
          (see the phiX_db.loc example file)
          
    c] Hyper variable in length amplicon databank
        - Upload databank from http://genoweb.toulouse.inra.fr/frogs_databanks/HVL
        - Extract databank.
        - To use this databank, you need to create a .loc file named
          'HVL_db.loc'. The path provided must be the '.fasta'.
          (see the HVL_db.loc example file)
          
### 9. Tools images
    The tools help contain images. These images must be in galaxy images
    static folder.
        ln -s <FROGS_PATH>/img <GALAXY_DIR>/static/images/tools/frogs


### 10. Use PEAR as reads merge software in preprocess
    PEAR is one of the most effective software for read pair merging, but as its licence is not free for private use, we can not distribute it in FROGS.
    If you work in an academic lab on a private Galaxy server, or if you have payed your licence you can use PEAR in FROGS preprocess.
    For that you need to:
    * have PEAR in your PATH or in the FROGS libexec directory
    * use --merge-software option for command line use
    * add PEAR in the FROGS preprocess Galaxy wrapper (<FROGS_DIR>/tools/preprocess/preprocess.xml): blok lines 117 tot 125 and block line 150 to 158 : 
    
    <conditional name="merge_software_type">
        <param name="merge_software" type="select" label="Merge software" help="Select the software to merge paired-end reads.">
            <option value="vsearch" selected="true">Vsearch</option>
            <option value="flash">Flash</option>
            <option value="pear">PEAR</option>
        </param>
        <when value="flash">
            <param name="expected_amplicon_size" type="integer" label="Expected amplicon size" help="Maximum amplicon length expected in approximately 90% of the amplicons." value="" />
        </when>
    </conditional>
        

## Troubleshooting
### Abnormal increase memory consumption with CPU number
    With certain old versions of glibc the virtual memory used by CPU is
    multiplicative.
    Nb CPUs   expected RAM consumtion   observed RAM consumption
       1               1Gb                       1Gb
       2               2Gb                     2*2Gb
       3               3Gb                     3*3Gb
       4               5Gb                     4*4Gb
    The parameters memory and CPU provided in examples take into account 
    this problem.

### Abnormal threads consumption in RDPClassifier
    With large database like silva NR the RDPClassifier opens automatically
    a large number of threads. These threads use all the available CPU
    ressources. This is not an acceptable behaviour in multi-user context.
    To prevent this behaviour the tool 'affiliation_OTU' uses taskset to
    force RDPClassifier to run only on the specified number of CPUs. The
    number of threads is not changed but the CPU consumption is controled.


## License
    GNU GPL v3


## Copyright
    2015 INRA


## Citation
    Please cite the **FROGS article**: *Escudie F., et al. Bioinformatics, 2018. FROGS: Find, Rapidly, OTUs with Galaxy Solution.*


## Contact
    frogs@inra.fr
=======
# Contact
frogs@inra.fr
>>>>>>> dev
