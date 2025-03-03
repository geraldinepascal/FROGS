#!/bin/bash
nb_cpu=$1
java_mem=$2
out_dir=$3

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
demultiplex.py --input-R1 data/demultiplex_test2_R1.fq.gz --input-R2 data/demultiplex_test2_R2.fq.gz --input-barcode data/demultiplex_barcode.txt \
    --mismatches 1 --end both \
    --output-demultiplexed $out_dir/demultiplexed.tar.gz --output-excluded $out_dir/undemultiplexed.tar.gz \
    --log-file $out_dir/demultiplex.log --summary $out_dir/demultiplex_summary.txt


if [ $? -ne 0 ]
then
	echo "Error in demultiplex " >&2
	exit 1;
fi

echo "Step denoising 16S vsearch `date`":


denoising.py illumina \
 --process swarm \
 --min-amplicon-size 44 --max-amplicon-size 490 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
 --R1-size 267 --R2-size 266 --merge-software vsearch \
 --nb-cpus $nb_cpu --mismatch-rate 0.15 \
 --input-archive data/test_dataset.tar.gz \
 --output-fasta $out_dir/01-denoising-swarm-vsearch.fasta \
 --output-biom $out_dir/01-denoising-swarm-vsearch.biom \
 --html $out_dir/01-denoising-swarm-vsearch.html \
 --log-file $out_dir/01-denoising-swarm-vsearch.log
 
 echo "Step denoising 16S vsearch swarm denoising and distance 3 `date`":


denoising.py illumina \
 --process swarm \
 --min-amplicon-size 44 --max-amplicon-size 490 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
 --R1-size 267 --R2-size 266 --merge-software vsearch \
 --nb-cpus $nb_cpu --mismatch-rate 0.15 \
 --input-archive data/test_dataset.tar.gz \
 --output-fasta $out_dir/01-denoising-swarm-dd3-vsearch.fasta \
 --output-biom $out_dir/01-denoising-swarm-dd3-vsearch.biom \
 --html $out_dir/01-denoising-swarm-dd3-vsearch.html \
 --log-file $out_dir/01-denoising-swarm-dd3--vsearch.log \
 --denoising --distance 3 

echo "Step denoising 16S pear `date`":

denoising.py illumina \
 --process swarm \
 --min-amplicon-size 44 --max-amplicon-size 490 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
 --R1-size 267 --R2-size 266 --merge-software pear \
 --nb-cpus $nb_cpu --mismatch-rate 0.15 \
 --input-archive data/test_dataset.tar.gz \
 --output-fasta $out_dir/01-denoising-swarm-pear.fasta \
 --output-biom $out_dir/01-denoising-swarm-pear.biom \
 --html $out_dir/01-denoising-swarm-pear.html \
 --log-file $out_dir/01-denoising-swarm-pear.log 

echo "Step denoising: dada2 keep-unmerged `date`"

denoising.py illumina  \
 --process dada2 --keep-unmerged \
 --input-archive data/verysmallITS.tar.gz \
 --min-amplicon-size 50 --max-amplicon-size 1000 --merge-software vsearch \
 --five-prim-primer TAGACTCGTCAHCGATGAAGAACGYRG --three-prim-primer GCATATCAATAAGCGSAGGAA \
 --R1-size 300 --R2-size 300  --nb-cpus $nb_cpu \
 --output-fasta $out_dir/01-denoising-dada2-clusters.fasta \
 --output-biom $out_dir/01-denoising-dada2-clusters.biom \
 --html $out_dir/01-denoising-dada2.html \
 --log-file $out_dir/01-denoising-dada2.log

echo "Step denoising: preprocess only `date`"

denoising.py illumina  \
 --process preprocess-only \
 --input-archive data/verysmallITS.tar.gz \
 --min-amplicon-size 50 --max-amplicon-size 1000 --merge-software vsearch \
 --five-prim-primer TAGACTCGTCAHCGATGAAGAACGYRG --three-prim-primer GCATATCAATAAGCGSAGGAA \
 --R1-size 300 --R2-size 300  --nb-cpus $nb_cpu \
 --output-fasta $out_dir/01-prepro-only-clusters.fasta \
 --output-biom $out_dir/01-prepro-only-clusters.biom \
 --html $out_dir/01-prepro-only.html \
 --log-file $out_dir/01-prepro-only.log

echo "Step remove_chimera `date`"

remove_chimera.py \
 --input-fasta $out_dir/01-denoising-swarm-vsearch.fasta \
 --input-biom $out_dir/01-denoising-swarm-vsearch.biom \
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
 --replicate_file data/replicates_file.tsv \
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
--impacted $out_dir/07-impacted_OTU_masked.tsv \
--impacted-multihit $out_dir/07-impacted_OTU_masked_multihit.tsv \
--log-file $out_dir/07-affiliation_filter_maskMode.log \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 150 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
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
--impacted $out_dir/07-impacted_OTU_deleted.tsv \
--impacted-multihit $out_dir/07-impacted_OTU_deleted_multihit.tsv \
--log-file $out_dir/07-affiliation_filter_delMode.log \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 150 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
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

echo "Step normalisation `date`"
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
    echo "Error in normalisation by min" >&2
    exit 1;
fi

# to reduce computing time for the others step
normalisation.py \
 --num-reads 100 \
 --input-biom $out_dir/08-affiliation_postprocessed.biom \
 --input-fasta $out_dir/08-affiliation_postprocessed.fasta \
 --output-biom $out_dir/09-normalisation.biom \
 --output-fasta $out_dir/09-normalisation.fasta \
 --html $out_dir/09-normalisation.html \
 --log-file $out_dir/09-normalisation.log
 
if [ $? -ne 0 ]
then
    echo "Error in normalisation by min" >&2
    exit 1;
fi

echo "Step cluster_stats `date`"

cluster_stats.py \
 --input-biom $out_dir/09-normalisation.biom \
 --html $out_dir/10-clustersStat.html \
 --log-file $out_dir/10-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stats" >&2
	exit 1;
fi

					
echo "Step affiliation_stats `date`"

affiliation_stats.py \
 --input-biom $out_dir/09-normalisation.biom \
 --html $out_dir/11-affiliationsStat.html \
 --log-file $out_dir/11-affiliationsStat.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species \
 --taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_stats" >&2
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
 --input-fasta $out_dir/09-normalisation.fasta \
 --input-biom $out_dir/09-normalisation.biom \
 --output-tree $out_dir/15-tree-mafft.nwk \
 --html $out_dir/15-tree-mafft.html \
 --log-file $out_dir/15-tree-mafft.log

if [ $? -ne 0 ]
then
	echo "Error in tree : mafft" >&2
	exit 1;
fi

echo "Step phyloseq_import_data `date`"

phyloseq_import_data.py  \
 --biomfile data/chaillou.biom \
 --samplefile data/sample_metadata.tsv \
 --treefile data/tree.nwk \
 --rdata $out_dir/16-phylo_import.Rdata \
 --html $out_dir/16-phylo_import.nb.html \
 --log-file $out_dir/16-phylo_import.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_import_data " >&2
	exit 1;
fi

echo "Step phyloseq_composition `date`"

phyloseq_composition.py  \
 --varExp EnvType --taxaRank1 Kingdom --taxaSet1 Bacteria --taxaRank2 Phylum --numberOfTaxa 9 \
 --rdata $out_dir/16-phylo_import.Rdata \
 --html $out_dir/17-phylo_composition.nb.html \
 --log-file $out_dir/17-phylo_composition.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_composition " >&2
	exit 1;
fi

echo "Step phyloseq_alpha_diversity `date`"

phyloseq_alpha_diversity.py  \
 --varExp EnvType \
 --rdata $out_dir/16-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --alpha-out $out_dir/18-phylo_alpha_div.tsv \
 --html $out_dir/18-phylo_alpha_div.nb.html \
 --log-file $out_dir/18-phylo_alpha_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_alpha_diversity " >&2
	exit 1;
fi

echo "Step phyloseq_beta_diversity `date`"

phyloseq_beta_diversity.py  \
 --varExp EnvType --distance-methods cc,unifrac \
 --rdata $out_dir/16-phylo_import.Rdata \
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
 --varExp EnvType --ordination-method MDS \
 --rdata $out_dir/16-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/20-phylo_structure.nb.html \
 --log-file $out_dir/20-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_structure " >&2
	exit 1;
fi

echo "Step phyloseq_clustering `date`"

phyloseq_clustering.py  \
 --varExp EnvType \
 --rdata $out_dir/16-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/21-phylo_clustering.nb.html \
 --log-file $out_dir/21-phylo_clustering.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_clustering " >&2
	exit 1;
fi

echo "Step phyloseq_manova `date`"

phyloseq_manova.py  \
 --varExp EnvType \
 --rdata $out_dir/16-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/22-phylo_manova.nb.html \
 --log-file $out_dir/22-phylo_manova.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_manova " >&2
	exit 1;
fi


echo "Step deseq2_preprocess `date`"
echo "DESeq2 asv abundances"

../app/deseq2_preprocess.py \
 --data $out_dir/16-phylo_import.Rdata \
 --analysis ASV \
 --log-file $out_dir/23-deseq2_preprocess_otu.log \
 --out-Rdata $out_dir/23-deseq2_preprocess_otu.Rdata \
 --var EnvType

echo "DESeq2 function abundances"

../app/deseq2_preprocess.py \
 --samplefile data/sample_metadata.tsv \
 --input-functions data/frogsfunc_functions_unstrat_EC.tsv \
 --analysis FUNCTION \
 --log-file $out_dir/23-deseq2_preprocess_func.log \
 --out-Rdata $out_dir/23-deseq2_preprocess_func.Rdata \
 --out-Phyloseq $out_dir/23-phyloseq_functions.Rdata \
 --var EnvType

if [ $? -ne 0 ]
then
	echo "Error in deseq2_preprocess " >&2
	exit 1;
fi


echo "Step deseq2_visualisation `date`"

echo "DESeq2 ASV abundances"
../app/deseq2_visualisation.py \
 --abundanceData $out_dir/16-phylo_import.Rdata \
 --analysis ASV \
 --dds $out_dir/23-deseq2_preprocess_otu.Rdata \
 --log-file $out_dir/24-deseq2_visualisation_otu.log \
 --html $out_dir/24-deseq2_visualisation_otu.nb.html \
 --var EnvType --mod1 BoeufHache --mod2 SaumonFume

echo "DESeq2 function abundances"
../app/deseq2_visualisation.py \
 --abundanceData $out_dir/23-phyloseq_functions.Rdata\
 --analysis FUNCTION \
 --dds $out_dir/23-deseq2_preprocess_func.Rdata \
 --log-file $out_dir/24-deseq2_visualisation_func.log \
 --html $out_dir/24-deseq2_visualisation_func.nb.html \
 --var EnvType --mod1 BoeufHache --mod2 SaumonFume

if [ $? -ne 0 ]
then
	echo "Error in deseq2_visualisation " >&2
	exit 1;
fi

echo "Completed with success"
