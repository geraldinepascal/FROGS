#Find Rapidly OTUs with Galaxy Solution

## Description
    FROGS is a galaxy/CLI workflow designed to produce an OTU count matrix 
    from high depth sequencing amplicon data.
    This workflow is focused on:
        - User-friendliness with the integration in galaxy and lots of rich 
          graphic outputs
        - Accuracy with a clustering without global similarity threshold, the
          management of multi-affiliations and management of separated PCRs
          in the chimera removal step
        - Speed with fast algorithms and an easy to use parallelisation
        - Scalability with algorithms designed to support the data growth


## Convenient input data
    Legend for the next schemas:
        .: Complete nucleic sequence
        !: Region of interest
        *: PCR primers
    
    Paired-end classical protocol:
        In the paired-end protocol R1 and R2 must share a nucleic region. 
        For example the amplicons on V3-V4 regions can have a length between
        350 and 500nt; with 2*300pb sequencing the overlap is between 250nt
        and 100nt. 
        From:                                    To:
         rDNA .........!!!!!!................    ......!!!!!!!!!!!!!!!!!!!.....
         Ampl      ****!!!!!!****                  ****!!!!!!!!!!!!!!!!!!!****
           R1      --------------                  --------------
           R2      --------------                               --------------
    
        The maximum overlap between R1 and R2 can be the complete overlap.
            Inconvenient maximum overlap:
            R1    --------------
            R2   --------------
        In this case it is necessary to trim R1 and R2 ends before the process.
    
        The minimum overlap between R1 and R2 can have 15nt. With less the 
        overlap can be incorrect.
              
    Single-end classical protocol:
        rDNA .........!!!!!!................
        Ampl      ****!!!!!!****
        Read      --------------

    Custom protocol
        rDNA .....!!!!!!!!!!!!!!............
        Ampl      ****!!!!!!****
        Read      --------------       

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
        Version: 1.1.3
        Named as: vsearch
        Tools: preprocess and remove_chimera
        Download: https://github.com/torognes/vsearch

    flash
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
                  
    vsearch
        Version: 1.1.3
        Named as: vsearch
        Tools: remove_chimera
        Download: https://github.com/torognes/vsearch
        
    swarm
        Version: 2.1.1
        Named as: swarm
        Tools: clustering
        Download: https://github.com/torognes/swarm

    ITSx
        Version : 1.0.11
        Named : ITSx
        Tools : preprocess
        Download : http://microbiology.se/software/itsx/
        Remark : ITSx_db folder need to be in the PATH or in <FROGS_PATH>/libexec
                 depend on HMMER 3 or later (only for hmmpress and hmmscan)
                 if ITSx test command line failed it's may be due to a difference in HMMER version used to prepare HMM models: hmmpress the ITSx_db/HMMs/*.hmm with your own HMMER version

    swarm
        Version: 2.1.1
        Named as: swarm
        Tools: clustering
        Download: https://github.com/torognes/swarm

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
        Version : = 3.4.0
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

    Dependencies must be accessible in <FROGS_PATH>/lib.

    Phyloseq-extended
        Version : v0.99
        Tools : all FROGSSTAT tools
        Installation : # https://github.com/mahendra-mariadassou/phyloseq-extended/releases
                       untar archive and copy content of folder "phyloseq-extended/" in <FROGS_PATH>/lib



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
        Step preprocess Fri Apr  8 11:44:10 CEST 2016
        Step clustering Fri Apr  8 11:45:32 CEST 2016
        Step remove_chimera Fri Apr  8 11:46:02 CEST 2016
        Step filters Fri Apr  8 11:47:11 CEST 2016
        Step ITSx Fri Apr  8 11:47:11 CEST 2016
        Step affiliation_OTU Fri Apr  8 11:47:12 CEST 2016
        Step clusters_stat Fri Apr  8 11:47:18 CEST 2016
        Step affiliations_stat Fri Apr  8 11:47:20 CEST 2016
        Step biom_to_tsv Fri Apr  8 11:47:40 CEST 2016
        Step biom_to_stdBiom Fri Apr  8 11:47:41 CEST 2016
        Step tsv_to_biom Fri Apr  8 11:47:42 CEST 2016
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

### 9. Tools images
    The tools help contain images. These images must be in galaxy images
    static folder.
        ln -s <FROGS_PATH>/img <GALAXY_DIR>/static/images/tools/frogs


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
    Escudie F., Auer L., Bernard M., Cauquil L., Vidal K., Maman S.,
    Mariadassou M., Hernadez-Raquet G., Pascal G., 2015. FROGS: Find Rapidly
    OTU with Galaxy Solution. In: The environmental genomic Conference, 
    Montpellier, France


## Contact
    frogs@toulouse.inra.fr
