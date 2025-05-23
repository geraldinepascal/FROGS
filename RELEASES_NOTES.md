# v5.0.2 [2025-05]
* Issue #77: remove max-len parameter for artificial combined sequences
* Upgrade mafft and blast versions

# v5.0.1 [2025-03]

### Bugs fixed
* Issue #75: unable to use denoising.py with swarm and --denoising option with --distance > 1
* splitbc.pl: fix R2 truncation
* denoising.py: Allow vsearch to deal with AVITI reads (higher max quality value allowed)

# v5.0.0 [2024-06]

### Tool added
* **FROGS denoising**: new tools replacing FROGS preprocess and FROGS clustering, allowing to perform the preprocess step and to chose between swarm and dada2 to produce ASVs.

### Tools removed
* **FROGS preprocess**: its functionality is integrated in the new FROGS denoising tool
* **FROGS clustering**: its functionality is integrated in the new FROGS denoising tool

### Global changes
* remove short version of parameters and standardization of parameter names
* add version and tool name in HTML reports
* add subject length information (slen) and related perc_subject_coverage in BIOM metadata for each affiliation
* remove and rename remaining tools/code with OTU term  
* possibility to swith color theme in HTML reports
  
### Functions added
* FROGS denoising:
  * add dada2 tool for Illumina paired-end data and PacBio long-reads data
  * add the possibility to perform a dereplication only instead of a clustering or denoising process
  * FROGS cluster stats tool integrated automatically in the HTML report
* FROGS demultiplex:
  * speed-up the way to check if file is gzipped
* FROGS taxonomic affiliation:
  * FROGS affiliation stats tool integrated automatically in the HTML report
  * Add subject length in BIOM
* FROGS remove chimera:
  * No more possibility to give a count matrix as input file
  * FROGS cluster stats tool integrated automatically in the HTML report
* FROGS itsx
  * FROGS cluster stats tool integrated automatically in the HTML report
* FROGS cluster stats
  * Hierarchical clustering no more done by default (option --hierarchical-clustering added)
* FROGS cluster filters
  * FROGS cluster stats tool integrated automatically in the HTML report
* FROGS phyloseq structure: use the same colors for samples in all ordination plots

### Bugs fixed
* Standardization of empty affiliations in BIOM files (empty list instead of None)
  * impacts FROGS affiliation filters, FROGS biom to tsv, FROGS tsv to biom, FROGS cluster filters, FROGS taxnomic affiliation, libexec/biomTools.py
* Phyloseq import and DESeq preprocess: allow files beginning with digits
* Phyloseq import: remove empty samples and gives the information in the HTML report
* DESeq preprocess: check samples names consistency between function abundance tsv file and sample_metadata
* DESeq visualization
  * remove empty samples and taxa before plotting the heatmap
  * correct heatmap where there is only one ASV or Function differentially expressed
  * correct ipath3 associated color (higher = darker) and reduce width difference
  * correct associated color of LFC>0 or LFC<0 in ipath3 images
* libexec/launch_hsp.py: correct path for tmp files
* ASV without affiliation is deleted if a threshold on a blast metric is set in affiliation_filters.py
* Affiliation_filters: Error message and exit if an ASV with a taxonomy has missing metrics (no data) because of a manual update
* FROGSfunc functions: correction of the location of temporary files 


# v4.1.0 [2023-03]

### Tools modified:
* FROGSFunc copynumber removed
* FROGSFunc placeseq : integrates FROGSFunc marker copynumber part
* FROGSFunc function : integrates FROGSFunc function copynumber part

### Function added
* Preprocess: add a longread sequencer option to deal with longreads
* Itsx : add mutual exclusion between --region and --check-its-only
* Affiliation_filter : add a –-keep-blast-taxa option : Taxon list to keep in Blast affiliations. All others affiliations will be masks/delete.
* DESeq preprocess : adapt the tool to deal with function abundance table (output from FROGSFunc Function) in TSV format
* DESeq visualisation : 
  * adapt the tool to functionnal abundance differential analysis
  * add Ipath3 visualisation functionnal abundance differential analysis
* FROGSFunc placeseq : 
  * add %identity and %coverage between amplicon sequences and Picrust reference sequences.
  * add a line chart of NSTI values versus blasts %identity and %coverage values.
* FROGSFunc function : 
  * in addition to NSTI, add filter on %identity and %coverage between amplicon sequences and Picrust reference sequences.
  * add a star plot to visualise the filter impact on taxonomies kept
  * parallelization of function abundance prediction on different function database (MetaCyc, KEGG, COG, ...), reducing the calculation time by more than half using all databases.

## Bug fixed
* Remove_chimera : deal with empty sample to compute proportion of sequence kept in the report. Proportion is now set to NA.
* Affiliation_filter : 
  * correctly find number of taxonomical rank
  * remove OTU without affiliation when using blast metrics filters

# v4.0.1 [2022-06]

### Bug fixed
* frogsfunc_placeseqs : 
  * repare html link in PICRUSt2 closest ID (JGI) column
  * add missing genomes in JGI_ID_to_taxonomy.txt file
  * deal with empty FROGS affiliation
  * add exception to avoid using sepp tools with ITS or 18S amplicon
* frogsfunc_function:
  * deal with function that are not associated with database link (picrust trait PHENO)
* Affiliation_filter: correctly find the number of taxonomical rank (in cas of empty affiliation in the first cluster)
* biom_to_stdBiom.py: deal with empty affiliations (add "Unclassified" for each taxonomic rank)

# v4.0.0 [2022-05]

### Tools added:

PICRUSt2 is a software for predicting functional abundances based only on marker gene sequences. This tool is integrated in S suite as FROGSFunc tools. They are splittedto 4 steps :
 * frogsfunc_placeseqs : places the OTUs into a reference phylogenetic tree.
 * frogsfunc_copynumbers : predicts marker and function copy number of each OTU.
 * frogsfunc_functions : calculates functions abundances in each sample.
 * frogsfunc_pathways :  calculates pathway abundances in each sample.

### Installation note:

As PICRUSt2 currently relies on a different R version, please install `frogsfunc-conda-requirements.yaml` and activate this environment before using FROGSFUNC tools. 

### Function added

  * Normalisation : 
    * add "Sampling by the number of sequences of the smallest sample" sampling method. This method automatically detects the sample with the smallest number of sequences, and samples all other samples with that number.  
    * If you chose "Select a number of reads" sampling method, you may or not activate "Remove samples" option. If it's activated, samples whose total number of sequences is lower than the specified number, will be removed from the abundance table. If the option is disabled, the samples will be kept in the analysis but with a number of sequences lower than the specified number (the total number of the sample). 
  * Otu_filter : add "Replicate identification" Minimum prevalence method. It allows to keep OTUs present in miniam replication proportions in at least one group (must be a proportion between 0 and 1). 
  * Affiliation_stat : add OTU rarefaction curves in HTML, in addition to the previously existing taxonomic ranks.
  * Remove_chimera : add "% Clusters kept" and "% Cluster abundance kept" in HTML chimera detection by sample table.

### Bug corrected:

* Affiliation_OTU : do not perform Needlall alignment if reduced reference database constructing by blasting R1 and R2 part of FROGS_combined OTU sequences is empty


# v3.2.3 [2021-05]

### Installation note:

When installing FROGS and its dependencies via conda, please refer de `frogs-conda-requirements.yaml` available on the master branch of the github repository ( https://github.com/geraldinepascal/FROGS/blob/master/frogs-conda-requirements.yaml ) instead of the one available in the release archive.
Indeed we may update dependencies versions without changing anything to the FROGS code so whithout making a new release.

### Bug fixed

* DESeq2 visualisation : correctly identify name of reference condition

# v3.2.2 [2021-04]

### ModificationsClustering

* Preprocess: use maxdiffpct instead of maxdiffs in vsearch fastq_mergepairs command line, and recommand 2.17.0 version. 
* DESeq2 : 
  * rename tool in DESeq2 visualisation (with s instead of z)
  * improve filter in datatable
  * change color
  * add padj threshold in MAplot
* ITSx : add organims model option (it was restrict to Fungi, take care of increase computing time
* OTU affiliation : sort blast affiliations in biom by taxonomy
* Clusters stat : add precision in HTML
* Remove chimera : add precision in HTML, and rename table columns names
* Affiliation Filter : add precision in HTML
* Various tools:
  * add taxonomic rank consistency between user declaration and input files (reference database, biom)
  * correct typo

### Bug fixed

* DESeq2 visualisation : 
  * add intermediates_dir argument in Rscript command
  * debug pie charts color attribution
* Normalisation : correct bug when calculating number of OTU by sample
* ITSx : correct stderr scanning
* Affiliation filter : correct bug in OTU filter by sample and by filter


# v3.2.1 [2021-02-22]

### Bug fixed
* correct reversecomp function by taking into account S and W IUPAC nucleotides

# v3.2 [2021-01-13]

### Tools added:
  * DESeq2 preprocess : Compute differential abundancy analysis
  * DESeq2 visualization : Create table and plots to explore and illustrate the differential abundant OTUs
  * Filters has been splitted into to new tools : FROGS OTU Filters and FROGS Affiliations Filters. 
    * FROGS OTU Filters filters OTU on presence/absence, abundances and contamination as Filters did. For contamination research, user may now use a personnal multifasta contaminant reference.
    * FROGS Affiliation Filters delete OTU or mask affiliation that do not respect affiliation metrics criteria, or affiliated to undesirable (partial) taxon.

### Function added
  * Affiliation_postprocess : taxon-ignore option added, to ignore some taxon like "unknown species" during the aggreagation process. Multiple taxon may be provided as well as partial taxon, like "sp."
  * Preprocess : now accept input sequence file as Fasta format (format automatically detected) for already contiged input.
  * Affiliations_stat : check that the number of rank name correspond to the number of ranks in input biom file
  * FROGSSTAT Phyloseq import : check sample names consistency between sample metadata and input biom file
  * FROGS OTU Filters filters OTU on presence/absence, abundances and contamination as Filters did.


### Bug fixed
  * Affiliation_postprocess : correctly compare de %coverage and the coverage threshold.
  * FROGSSTAT Structure : plot_heatmap now take into account ordination method and dissimilarity matrix
  * addAffiliation2biom : do not split partial description from blast reference ID
  * tsv_2_biom now keep initial OTU order (Cluster_1 is the most abundant one and Cluster_X the less abundant one)

### Other improvements
  * All python scripts are now in python 3
  * FROGSSTAT Phyloseq and FROGSSTAT DESeq now generate notebook_html instead of classical HTML output file. This facilitates code maintenance.

# v3.1 [2018-01-08]

### Bug fixed
  * Tsv_to_Biom : manage quotes added by Excel
  * Affiliation_postprocess : bug when references are not in the troncated amplicon database.

### Modifications:
  * Tree do no longer support Pynast alignment thanks to a template file
  * XML wrappers, images and loc files are externalized to [FROGS-wrappers](https://github.com/geraldinepascal/FROGS-wrappers) (FROGS is now available on Toolshed and conda )

# v3.0.0 [2018-10-10]
### Tools added:
  * ITSx : tool available for selecting and trimming ITS sequences based on ITSx tool
  * Affiliation Postprocess : resolve ambiguities due to inclusiv ITS, and aggregated OTU based on 
    taxonomic affiliations

### Functions added:
  * Preprocess : keep and filter non overlapped reads (particularly important for amplicon polymorphe in length)
  * Affiliation : non overlapped read generate FROGS_combined OTU that are affiliated thanks to Needl (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html) global alignment

### Modifications:
  * General : Use GALAXY_SLOTS variable in XML. 
              Upadate version dependancies
              Remove redirection to /dev/null
              Update test/test.sh
              Update XML doc
              Update graphical charter of HTML outputs
              HTML reports now work on secured server by using of https instead of http to all JavaScript and CSS URLs (https://github.com/geraldinepascal/FROGS/issues/34).
  * Preprocess : add VSEARCH and PEAR (https://sco.h-its.org/exelixis/web/software/pear/ only for command line usage because of PEAR licence) options to merge reads
  * Clustering : deal with artificially combined reads
  * FROGSSTAT tools do not call anymore external script from urls. Additionnal library dependancies 
    added : https://github.com/mahendra-mariadassou/phyloseq-extended/releases
  * FROGSSTAT import : can take FROGS specific biom file
  * FROGS TSV_to_Biom : keeps original type of metadata (float, int or list instead of string)
                        remove quotes added by Excel when exporting abundance tsv file in text file

### Bugs fixes:
  * Phyloseq import: bug when Tree is None
  * libexec/derepSamples : bug with temporary files name

# v2.0.0  [2017-08-08]
### Tools added : 
  * Tree : perform phylogenetic tree reconstruction based on Mafft or Pynast follow by Fasttree and phangorn
  * FROGS Phyloseq Import Data : import data from 3 files: biomfile, samplefile, treefile into a phyloseq R object.
  * FROGS Phyloseq Composition Visualization with bar plot and composition plot
  * FROGS Phyloseq Alpha Diversity with richness plot
  * FROGS Phyloseq Beta Diversity distance matrix
  * FROGS Phyloseq Structure
  * FROGS Phyloseq Clustering
  * FROGS Phyloseq Manova 

### libexec program added:
  * rooted_tree.R : Rscript to root FastTree tree. (used by Tree)

### Bugs fixes:
  * Preprocess : min overlap at least equal to 1
  * biom2tsv : not working with stdBiom containing RDP affiliation, not working when emtpy metadata

### Functions added:
  * Preprocess: add Flash mismatch rate option

# v1.4.0  [2017-02-04]
### Bugs fixes:
  * Preprocess: error in final dereplication with hudge number of samples.
  * Remove_chimera: error when using library Queue and hudge number of samples.
  * Clusters_stat: error with empty samples in hierarchical clustering.
  * Filters: error when only the filter on contamination is used.
  * Filters: bug when using other filters than abundance (check parameter when None).
  * Tsv2Biom : bug fix when using a tsv file comming from a standard biom file
  * Affiliations_stat : bug in rarefaction step computation when sample are empty

### Functions added:
  * Preprocess: new amplicon length graph.
  * Clustering: reduce memory consumption and execution time for the step swarm2biom.
  * Affiliations_stat: more details in alignment heatmap.


# v1.3.0  [2016-04-18]
### Update notes:
This section contain actions to execute if you update a previous FROGS version to this version.

  * The directory lib has been created and librairies have been renamed. In installation step 2 if you have added `FROGS/bin` in the PYTHONPATH you can remove this path and add `FROGS/lib` instead.
  * The directory libexec has been created. In installation step 2 if you have added `FROGS/bin` in the PATH you can remove this path and add `FROGS/libexec` instead.
  * The tools are moved in app directory. If you use FROGS in galaxy you must change the tool_conf.xml:

        ...
        <section id="FROGS_wrappers" name="FROGS">
            <tool file="FROGS/app/upload_tar.xml" />
            <tool file="FROGS/app/demultiplex.xml" />
            <tool file="FROGS/app/preprocess.xml" />
            <tool file="FROGS/app/clustering.xml" />
            <tool file="FROGS/app/remove_chimera.xml" />  
            <tool file="FROGS/app/filters.xml" />
            <tool file="FROGS/app/affiliation_OTU.xml" />
            <tool file="FROGS/app/clusters_stat.xml" />
            <tool file="FROGS/app/affiliations_stat.xml" />
            <tool file="FROGS/app/biom_to_stdBiom.xml" />
            <tool file="FROGS/app/biom_to_tsv/biom_to_tsv.xml" />
            <tool file="FROGS/app/tsv_to_biom.xml" />
            <tool file="FROGS/app/normalisation.xml" />
        </section>
        ...

### Bugs fixes:
  * Preprocess: missing values in lengths distribution graph.


# v1.2.0  [2016-04-06]
### Update notes:
This section contain actions to execute if you update a previous FROGS version to this version.

  * The repository structure has been changed: tools have been moved in sub-directories. It is necessary to re-do the following installation steps: '6. Add tools in galaxy' and '7. Set memory and parallelisation settings'.
  * The tool FROGS\_tsv\_to\_biom has been added. It is necessary to add the following line in galaxy tool_conf.xml: `<tool file="FROGS/tools/tsv_to_biom/tsv_to_biom.xml" />`.

### Bugs fixes:
  * Too large number of temp files in clustering tool.
  * Fix bug with BIOM without observation metadata in filters.
  * Fix bug with heatmap hover in affiliations_stat.

### Functions added:
  * Add tool to convert a TSV file in BIOM file.
  * Default BIOM becomes mono-line.
  * Reduce the BIOM file size (identation 4 -> 2 in pretty_print).
  * Change the repository structure (one directory by tool).
  * Use biomTools in normalisation.
  * Number of sampled sequences becomes required in normalisation tool.
  * Add better wrapping for sample name fields in preprocess.
  * Add check on minimum amplicon size value in preprocess.
  * Add vsearch version in remove chimera log.
  * Add cutadapt version in preprocess log.


# v1.1.0  [2015-11-30]
### Bugs fixes:
  * Fix bug with empty blast_taxonomy in BIOM.
  * Fix bug with checkall in chrome.
  * Fix bug with CSV extraction in HTML outputs.
  * Fix bug 'Exception in thread QueueFeederThread' in slow systems.

### Functions added:
  * Add '% kept' in preprocess.
  * Add checks on samples names in preprocess.
  * Add DenseData management in BIOM.
  * Add the biom1 datatype in tools.


# v1.0.0  [2015-09-18]
  First package.
