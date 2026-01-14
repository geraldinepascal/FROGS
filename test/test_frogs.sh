#!/bin/bash
nb_cpu=$1
java_mem=$2
out_dir=$3

# conda activate frogs@XX
# export PATH=~/workspace/FROGS_dev/app:$PATH

# Check parameters
if [ "$#" -ne 3 ]; then
    echo "ERROR: Illegal number of parameters." ;
    echo 'Command usage: test_frogs.sh <NB_CPU> <JAVA_MEM> <OUT_FOLDER>' ;
    exit 1 ;
fi

# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi

echo "Step demultiplex `date`"

demultiplex.py \
	--input-R1 data/formation_small.fastq.gz \
	--input-barcode data/formation.barcode.txt \
	--mismatches 1 --end both \
	--output-demultiplexed $out_dir/demultiplexed.tar.gz \
	--output-undemultiplexed $out_dir/undemultiplexed.tar.gz \
	--log-file $out_dir/demultiplex.log \
	--summary $out_dir/demultiplex_summary.txt 

if [ $? -ne 0 ]
then
	echo "Error in demultiplex " >&2
	exit 1;
fi

echo "Step reads_processing 16S vsearch swarm fastidious `date`":

reads_processing.py illumina \
	--process swarm --fastidious \
	--min-amplicon-size 44 --max-amplicon-size 490 \
	--five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
	--R1-size 267 --R2-size 266 --merge-software vsearch \
	--nb-cpus $nb_cpu --mismatch-rate 0.15 \
	--input-archive data/test_dataset.tar.gz \
	--output-fasta $out_dir/01-reads_processing-swarm-vsearch.fasta \
	--output-biom $out_dir/01-reads_processing-swarm-vsearch.biom \
	--output-compo $out_dir/01-reads_processing-swarm-vsearch_compo.tsv \
	--html $out_dir/01-reads_processing-swarm-vsearch.html \
	--log-file $out_dir/01-reads_processing-swarm-vsearch.log

if [ $? -ne 0 ]
then
	echo "Error in reads_processing 16S vsearch " >&2
	exit 1;
fi
 
 echo "Step reads_processing 16S vsearch swarm reads_processing and distance 3 `date`":

reads_processing.py illumina \
	--process swarm \
	--min-amplicon-size 44 --max-amplicon-size 490 \
	--five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
	--R1-size 267 --R2-size 266 --merge-software vsearch \
	--nb-cpus $nb_cpu --mismatch-rate 0.15 \
	--input-archive data/test_dataset.tar.gz \
	--output-fasta $out_dir/01-reads_processing-swarm-dd3-vsearch.fasta \
	--output-biom $out_dir/01-reads_processing-swarm-dd3-vsearch.biom \
	--output-compo $out_dir/01-reads_processing-swarm-dd3-vsearch_compo.tsv \
	--html $out_dir/01-reads_processing-swarm-dd3-vsearch.html \
	--log-file $out_dir/01-reads_processing-swarm-dd3--vsearch.log \
	--pre-clustering --distance 3 

if [ $? -ne 0 ]
then
	echo "Error in reads_processing 16S vsearch swarm reads_processing and distance 3 " >&2
	exit 1;
fi

echo "Step reads_processing: dada2 keep-unmerged `date`"

reads_processing.py illumina  \
	--process dada2 --keep-unmerged \
	--input-archive data/verysmallITS.tar.gz \
	--min-amplicon-size 50 --max-amplicon-size 1000 --merge-software vsearch \
	--five-prim-primer TAGACTCGTCAHCGATGAAGAACGYRG --three-prim-primer GCATATCAATAAGCGSAGGAA \
	--R1-size 300 --R2-size 300  --nb-cpus $nb_cpu \
	--output-fasta $out_dir/01-reads_processing-dada2-clusters.fasta \
	--output-biom $out_dir/01-reads_processing-dada2-clusters.biom \
	--html $out_dir/01-reads_processing-dada2.html \
	--log-file $out_dir/01-reads_processing-dada2.log

if [ $? -ne 0 ]
then
	echo "Error in reads_processing: dada2 keep-unmerged " >&2
	exit 1;
fi

echo "Step reads_processing: preprocess only `date`"

reads_processing.py illumina  \
	--process preprocess-only \
	--input-archive data/verysmallITS.tar.gz \
	--min-amplicon-size 50 --max-amplicon-size 1000 --merge-software vsearch \
	--five-prim-primer TAGACTCGTCAHCGATGAAGAACGYRG --three-prim-primer GCATATCAATAAGCGSAGGAA \
	--R1-size 300 --R2-size 300  --nb-cpus $nb_cpu \
	--output-fasta $out_dir/01-prepro-only-clusters.fasta \
	--output-biom $out_dir/01-prepro-only-clusters.biom \
	--html $out_dir/01-prepro-only.html \
	--log-file $out_dir/01-prepro-only.log

if [ $? -ne 0 ]
then
	echo "Error in reads_processing: preprocess only " >&2
	exit 1;
fi

echo "Step remove_chimera `date`"

remove_chimera.py \
	--input-fasta $out_dir/01-reads_processing-swarm-vsearch.fasta \
	--input-biom $out_dir/01-reads_processing-swarm-vsearch.biom \
	--output-fasta $out_dir/03-chimera.fasta \
	--output-biom $out_dir/03-chimera.biom \
	--html $out_dir/03-chimera.html \
	--log-file $out_dir/03-chimera.log \
	--nb-cpus $nb_cpu
 
if [ $? -ne 0 ]
then
	echo "Error in remove_chimera" >&2
	exit 1;
fi

echo "Step cluster_filters `date`"

cluster_filters.py \
	--min-abundance 0.00005 \
	--min-sample-presence 3 \
	--contaminant data/phi.fa \
	--nb-cpus $nb_cpu \
	--input-biom $out_dir/03-chimera.biom \
	--input-fasta $out_dir/03-chimera.fasta \
	--replicate-tsv data/replicates_file.tsv \
	--min-replicate-presence 0.5 \
	--output-fasta $out_dir/04-filters.fasta \
	--output-biom $out_dir/04-filters.biom \
	--excluded $out_dir/04-filters.excluded \
	--html $out_dir/04-filters.html \
	--log-file $out_dir/04-filters.log 

if [ $? -ne 0 ]
then
	echo "Error in cluster_filters" >&2
	exit 1;
fi

echo "Step itsx `date`"

itsx.py \
	--input-fasta $out_dir/04-filters.fasta \
	--input-biom $out_dir/04-filters.biom \
	--region ITS1 --nb-cpus $nb_cpu \
	--output-biom $out_dir/05-itsx.biom \
	--html $out_dir/05-itsx.html \
	--log-file $out_dir/05-itsx.log \
	--output-fasta $out_dir/05-itsx.fasta \
	--output-removed-sequences $out_dir/05-itsx-excluded.fasta

if [ $? -ne 0 ]
then
	echo "Error in ITSx" >&2
	exit 1;
fi

echo "Step taxonomic_affiliation `date`"

taxonomic_affiliation.py \
	--reference data/ITS1.rdp.fasta \
	--input-fasta $out_dir/04-filters.fasta \
	--input-biom $out_dir/04-filters.biom \
	--output-biom $out_dir/06-affiliation.biom \
	--html $out_dir/06-affiliation.html \
	--log-file $out_dir/06-affiliation.log \
	--nb-cpus $nb_cpu --java-mem $java_mem \
	--rdp

if [ $? -ne 0 ]
then
	echo "Error in taxonomic_affiliation" >&2
	exit 1;
fi

echo "Step affiliation_filters: masking mode `date`"

affiliation_filters.py \
	--input-biom $out_dir/06-affiliation.biom \
	--input-fasta $out_dir/04-filters.fasta \
	--output-biom $out_dir/07-affiliation_masked.biom \
	--html $out_dir/07-affiliation_masked.html \
	--impacted $out_dir/07-impacted_ASV_masked.tsv \
	--impacted-multihit $out_dir/07-impacted_ASV_masked_multihit.tsv \
	--log-file $out_dir/07-affiliation_filter_maskMode.log \
	--min-blast-length 150 \
	--min-blast-identity 100 \
	--min-blast-coverage 100 \
	--max-blast-evalue 1e-150 \
	--ignore-blast-taxa "g__Sarcodon" "s__Trichoderma" \
	--mask \
	--taxonomic-ranks Domain Phylum Class Order Family Genus Species  # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_filters: masking mode" >&2
	exit 1;
fi

echo "Step affiliation_filters: deleted mode `date`"

affiliation_filters.py \
	--input-biom $out_dir/06-affiliation.biom \
	--input-fasta $out_dir/04-filters.fasta \
	--output-biom $out_dir/07-affiliation_deleted.biom \
	--output-fasta $out_dir/07-affiliation_deleted.fasta \
	--html $out_dir/07-affiliation_deleted.html \
	--impacted $out_dir/07-impacted_ASV_deleted.tsv \
	--impacted-multihit $out_dir/07-impacted_ASV_deleted_multihit.tsv \
	--log-file $out_dir/07-affiliation_filter_delMode.log \
	--min-blast-length 150 \
	--min-blast-identity 100 \
	--min-blast-coverage 100 \
	--max-blast-evalue 1e-150 \
	--ignore-blast-taxa "g__Sarcodon" "s__Trichoderma" \
	--delete \
	--taxonomic-ranks Domain Phylum Class Order Family Genus Species  # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_filters: deleted mode" >&2
	exit 1;
fi

echo "Step affiliation_postprocess `date`"

affiliation_postprocess.py \
	--input-biom $out_dir/06-affiliation.biom \
	--input-fasta $out_dir/04-filters.fasta \
	--reference data/Unite_extract_ITS1.fasta \
	--output-biom $out_dir/08-affiliation_postprocessed.biom \
	--output-compo $out_dir/08-affiliation_postprocessed.compo.tsv \
	--output-fasta $out_dir/08-affiliation_postprocessed.fasta \
	--log-file $out_dir/08-affiliation_postprocessed.log

if [ $? -ne 0 ]
then
	echo "Error in affiliation_postprocess" >&2
	exit 1;
fi

echo "Step normalisation fixe num-reads delete sample `date`"

normalisation.py \
	--num-reads 25000 \
	--delete-samples \
	--input-biom $out_dir/08-affiliation_postprocessed.biom \
	--input-fasta $out_dir/08-affiliation_postprocessed.fasta \
	--output-biom $out_dir/09-normalisation_25K_delS.biom \
	--output-fasta $out_dir/09-normalisation_25K_delS.fasta \
	--html $out_dir/09-normalisation_25K_delS.html \
	--log-file $out_dir/09-normalisation_25K_delS.log
 
if [ $? -ne 0 ]
then
	echo "Error in normalisation 25K_delS" >&2
	exit 1;
fi

echo "Step normalisation fixe min_sample `date`"

normalisation.py \
	--sampling-by-min \
	--input-biom $out_dir/08-affiliation_postprocessed.biom \
	--input-fasta $out_dir/08-affiliation_postprocessed.fasta \
	--output-biom $out_dir/09-normalisation_by_min.biom \
	--output-fasta $out_dir/09-normalisation_by_min.fasta \
	--html $out_dir/09-normalisation_by_min.html \
	--log-file $out_dir/09-normalisation_by_min.log
 
if [ $? -ne 0 ]
then
    echo "Error in normalisation by min_sample" >&2
    exit 1;
fi

echo "Step normalisation fixe small num-reads`date`"

normalisation.py \
	--num-reads 100 \
	--delete-samples \
	--input-biom $out_dir/08-affiliation_postprocessed.biom \
	--input-fasta $out_dir/08-affiliation_postprocessed.fasta \
	--output-biom $out_dir/09-normalisation.biom \
	--output-fasta $out_dir/09-normalisation.fasta \
	--html $out_dir/09-normalisation.html \
	--log-file $out_dir/09-normalisation.log
 
if [ $? -ne 0 ]
then
    echo "Error in normalisation by fixe small num-reads" >&2
    exit 1;
fi

echo "Step cluster_asv_report `date`"

cluster_asv_report.py \
 --input-biom $out_dir/09-normalisation.biom \
 --html $out_dir/10-clusters-asv-report.html \
 --log-file $out_dir/10-clusters-asv-report.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stats" >&2
	exit 1;
fi
	
echo "Step affiliation_report `date`"

affiliation_report.py \
 --input-biom $out_dir/09-normalisation.biom \
 --html $out_dir/11-affiliation_report.html \
 --log-file $out_dir/11-affiliation_report.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species \
 --taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_report" >&2
	exit 1;
fi

echo "Step biom_to_tsv `date`"

biom_to_tsv.py \
 --input-biom $out_dir/09-normalisation.biom \
 --input-fasta $out_dir/09-normalisation.fasta \
 --output-tsv $out_dir/12-biom2tsv.tsv \
 --output-multi-affi $out_dir/12-biom2tsv-affiliation_multihit.tsv \
 --log-file $out_dir/12-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi

echo "Step biom_to_stdBiom `date`"

biom_to_stdBiom.py \
 --input-biom $out_dir/09-normalisation.biom \
 --output-biom $out_dir/13-affiliation_std.biom \
 --output-metadata $out_dir/13-affiliation_multihit.tsv \
 --log-file $out_dir/13-biom2stdbiom.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_stdBiom" >&2
	exit 1;
fi

echo "Step tsv_to_biom `date`"

tsv_to_biom.py \
 --input-tsv $out_dir/12-biom2tsv.tsv \
 --input-multi-affi $out_dir/12-biom2tsv-affiliation_multihit.tsv \
 --output-biom $out_dir/14-tsv2biom.biom \
 --output-fasta $out_dir/14-tsv2biom.fasta \
 --log-file $out_dir/14-tsv2biom.log 

if [ $? -ne 0 ]
then
	echo "Error in tsv_to_biom" >&2
	exit 1;
fi

echo "Step tree `date`"

tree.py \
 --nb-cpus $nb_cpu \
 --input-fasta $out_dir/04-filters.fasta \
 --input-biom $out_dir/06-affiliation.biom \
 --output-tree $out_dir/15-tree-mafft.nwk \
 --html $out_dir/15-tree-mafft.html \
 --log-file $out_dir/15-tree-mafft.log

if [ $? -ne 0 ]
then
	echo "Error in tree " >&2
	exit 1;
fi

echo "Step phyloseq_import_data `date`"

phyloseq_import_data.py  \
 --input-biom data/chaillou.biom \
 --sample-metadata-tsv data/sample_metadata.tsv \
 --tree-nwk data/tree.nwk \
 --output-phyloseq-rdata $out_dir/16-phylo_import.Rdata \
 --html $out_dir/16-phylo_import.nb.html \
 --log-file $out_dir/16-phylo_import.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_import_data " >&2
	exit 1;
fi

echo "Step phyloseq_composition `date`"

phyloseq_composition.py  \
 --var-exp EnvType --taxa-rank-1 Kingdom --taxa-set-1 Bacteria --taxa-rank-2 Phylum --number-of-taxa 9 \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata \
 --html $out_dir/17-phylo_composition.nb.html \
 --log-file $out_dir/17-phylo_composition.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_composition " >&2
	exit 1;
fi

echo "Step phyloseq_alpha_diversity `date`"

phyloseq_alpha_diversity.py  \
 --var-exp EnvType \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --output-alpha-tsv $out_dir/18-phylo_alpha_div.tsv \
 --html $out_dir/18-phylo_alpha_div.nb.html \
 --log-file $out_dir/18-phylo_alpha_div.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_alpha_diversity " >&2
	exit 1;
fi

echo "Step phyloseq_beta_diversity `date`"

phyloseq_beta_diversity.py  \
 --var-exp EnvType --beta-distance-methods cc unifrac \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata \
 --matrix-outdir $out_dir \
 --html $out_dir/19-phylo_beta_div.nb.html \
 --log-file $out_dir/19-phylo_beta_div.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_beta_diversity " >&2
	exit 1;
fi

echo "Step phyloseq_structure `date`"

phyloseq_structure.py  \
 --var-exp EnvType --ordination-method MDS \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata --beta-distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/20-phylo_structure.nb.html \
 --log-file $out_dir/20-phylo_structure.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_structure " >&2
	exit 1;
fi

echo "Step phyloseq_clustering `date`"

phyloseq_clustering.py  \
 --var-exp EnvType \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata --beta-distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/21-phylo_clustering.nb.html \
 --log-file $out_dir/21-phylo_clustering.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_clustering " >&2
	exit 1;
fi

echo "Step phyloseq_manova `date`"

phyloseq_manova.py  \
 --var-exp EnvType \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata --beta-distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/22-phylo_manova.nb.html \
 --log-file $out_dir/22-phylo_manova.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_manova " >&2
	exit 1;
fi

echo "Step deseq2_preprocess ASV `date`"

deseq2_preprocess.py \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata \
 --analysis-type ASV \
 --log-file $out_dir/23-deseq2_preprocess_asv.log \
 --output-deseq-rdata $out_dir/23-deseq2_preprocess_asv.Rdata \
 --var-exp EnvType

if [ $? -ne 0 ]
then
	echo "Error in deseq2_preprocess ASV " >&2
	exit 1;
fi

echo "Step deseq2_preprocess Function `date`"

deseq2_preprocess.py \
 --sample-metadata-tsv data/sample_metadata.tsv \
 --input-functions-abund data/frogsfunc_functions_unstrat_EC.tsv \
 --analysis FUNCTION \
 --log-file $out_dir/23-deseq2_preprocess_func.log \
 --output-deseq-rdata $out_dir/23-deseq2_preprocess_func.Rdata \
 --output-phyloseq-rdata $out_dir/23-phyloseq_functions.Rdata \
 --var-exp EnvType

if [ $? -ne 0 ]
then
	echo "Error in deseq2_preprocess Function " >&2
	exit 1;
fi

echo "Step deseq2_visualisation ASV `date`"

deseq2_visualisation.py \
 --phyloseq-rdata $out_dir/16-phylo_import.Rdata \
 --analysis-type ASV \
 --deseq-rdata $out_dir/23-deseq2_preprocess_asv.Rdata \
 --log-file $out_dir/24-deseq2_visualisation_asv.log \
 --html $out_dir/24-deseq2_visualisation_asv.nb.html \
 --var EnvType --mod1 BoeufHache --mod2 SaumonFume

if [ $? -ne 0 ]
then
	echo "Error in deseq2_visualisation ASV" >&2
	exit 1;
fi

echo "Step deseq2_visualisation Function `date`"

deseq2_visualisation.py \
 --phyloseq-rdata $out_dir/23-phyloseq_functions.Rdata\
 --analysis FUNCTION \
 --deseq-rdata $out_dir/23-deseq2_preprocess_func.Rdata \
 --log-file $out_dir/24-deseq2_visualisation_func.log \
 --html $out_dir/24-deseq2_visualisation_func.nb.html \
 --var EnvType --mod1 BoeufHache --mod2 SaumonFume

if [ $? -ne 0 ]
then
	echo "Error in deseq2_visualisation Function" >&2
	exit 1;
fi
echo "Completed with success"
