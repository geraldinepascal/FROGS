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
### 1. Install dependencies
	python interpreter
		Version: 2.7
		Tools: all
	python SciPy
		Tools: clusters_stat
	perl interpreter
		Version: 5
		Tools: demultiplex
	flash
		Version: >=1.2.8
		Named as: flash
		Tools: preprocess
		Download: http://sourceforge.net/projects/flashpage/files
	cutadapt
		Version: >=1.7
		Named as: cutadapt
		Tools: preprocess
		Download: https://github.com/marcelm/cutadapt
				  OR
				  https://pypi.python.org/pypi/cutadapt
	swarm
		Version: >=2.1.1
		Named as: swarm
		Tools: clustering
		Download: https://github.com/torognes/swarm
	vsearch
		Version: >=1.1.3
		Named as: vsearch
		Tools: remove_chimera
		Download: https://github.com/torognes/vsearch
	NCBI Blast+ blastn
		Version: >=2.2.29+
		Named as: blastn
		Tools: affiliation_OTU and filters
		Download: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
	RDPClassifier
		Version: -
		Name as: classifier.jar
		Tools: affiliation_OTU
		Download: https://github.com/rdpstaff/RDPTools
	taskset
		Version: util-linux-ng 2.17.2
		Name as: taskset
		Tools: affiliation_OTU
		Install: sudo apt-get install util-linux
				 OR
				 sudo yum install util-linux
    
### 2. Bin directory
	Softwares in 'bin' folder and above dependencies must be in 
	PATH and PYTHONPATH or in a 'bin' folder in the parent folder of the
	tools folder.
	  Example with FROGS binaries in bin folder:
		<FROGS_PATH>/
			tools/
				preprocess.py
				preprocess.xml
				...
			bin/
				biom2tsv.py
				biom.py
				...

### 3. Check intallation
	To check your installation you can type:
		cd <FROGS_PATH>/test
		bash test.sh <FROGS_PATH> <NB_CPU> <JAVA_MEM> <OUT_FOLDER>
	
	This test executes the FROGS tools in command line mode.
	Note:
		JAVA_MEM must be at least 4 (= 4Gb of RAM).
	Example:
		[user@computer:/home/user]$cd /home/user/frogs_git/test
		[user@computer:/home/user/frogs_git/test]$bash test.sh /home/user/frogs_git/ 2 4 results
		Step preprocess Thu Sep 24 08:34:01 CEST 2015
		Step clustering Thu Sep 24 08:35:38 CEST 2015
		Step remove_chimera Thu Sep 24 08:37:32 CEST 2015
		Step filters Thu Sep 24 08:38:27 CEST 2015
		Step affiliation_OTU Thu Sep 24 08:38:30 CEST 2015
		Step clusters_stat Thu Sep 24 08:38:44 CEST 2015
		Step affiliations_stat Thu Sep 24 08:38:51 CEST 2015
		Step biom_to_tsv Thu Sep 24 08:39:13 CEST 2015
		Completed with success

### 4. Add tools in galaxy
	Add the tools in galaxy tool_conf.xml.
	Example:
		...
		<section id="FROGS_wrappers" name="FROGS">
			<tool file="FROGS/tools/upload_tar.xml" />
			<tool file="FROGS/tools/demultiplex.xml" />
			<tool file="FROGS/tools/preprocess.xml" />
			<tool file="FROGS/tools/clustering.xml" />
			<tool file="FROGS/tools/remove_chimera.xml" />  
			<tool file="FROGS/tools/filters.xml" />
			<tool file="FROGS/tools/affiliation_OTU.xml" />
			<tool file="FROGS/tools/biom_to_tsv.xml" />
			<tool file="FROGS/tools/clusters_stat.xml" />
			<tool file="FROGS/tools/affiliations_stat.xml" />
			<tool file="FROGS/tools/biom_to_stdBiom.xml" />
			<tool file="FROGS/tools/normalisation.xml" />
		</section>
		...

### 5. Set memory and parallelisation settings
	If you have more than one CPU, it is recommended to increase the number
	of CPUs used by tools.
	All the CPUs must be on the same computer/node.

	a] Specifications  
		Tool            RAM/CPU     Minimal RAM     Configuration example
		affiliation	    -              20 Gb         30 CPUs and 300 GB
		chimera         3 Gb           5 Gb          12 CPUs and 36 GB
		clustering	    -              10 Gb         16 CPUs and 60 GB
		preprocess	    8 Gb            -            12 CPUs and 96 GB

	b] Change the number of CPUs used
		Each tool with parallelisation possibilities contains in its XML an
		hidden parameter to set the number of used CPUs.
		Example for 16 CPUs:
		   <param name="nb_cpus" type="hidden" label="CPU number" help="The maximum number of CPUs used." value="1" />
		   Is changed to :
		   <param name="nb_cpus" type="hidden" label="CPU number" help="The maximum number of CPUs used." value="16" />

	c] Change the tool launcher configuration
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
				<tool id="FROGS_affiliation_OTU" destination="FROGS_affiliation_OTU_job"/>
			</tools>

### 6. Upload and configure the databanks
	a] Assignation databank
		- Upload databanks and indexes from http://genoweb.toulouse.inra.fr/frogs_databanks/assignation
		- Extract databanks.
		- To use these databank, you need to create a .loc file named
		  'frogs_db.loc'. The path provided must be the '.fasta'.
	b] Contaminant databank
		- Upload databank and indexes from http://genoweb.toulouse.inra.fr/frogs_databanks/contaminants
		- Extract databank.
		- To use this databank, you need to create a .loc file named
		  'phiX_db.loc'. The path provided must be the '.fasta'.

### 7. Tools images
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

## 

## Citation
    Escudie F., Auer L., Cauquil L., Vidal K., Maman S., Mariadassou M., 
    Hernadez-Raquet G., Pascal G., 2015. FROGS: Find Rapidly OTU with Galaxy 
    Solution. In: (Ed), The JOBIM 2015 Conference, July 6th to 9th, 
    Clermont-Ferrand, France


## Contact
    frogs@toulouse.inra.fr
