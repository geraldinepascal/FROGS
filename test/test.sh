#!/bin/bash
frogs_dir=$1
nb_cpu=$2
java_mem=$3


# Set ENV
export PATH=$frogs_dir/bin:$PATH
export PYTHONPATH=$frogs_dir/bin:$PYTHONPATH


# Clean previous results
if [ ! -d "$frogs_dir/test/results" ]
then
    mkdir $frogs_dir/test/results
fi
rm -f $frogs_dir/test/results/*


echo "Step preprocess `date`"

$frogs_dir/tools/preprocess.py illumina \
 --min-amplicon-size 380 --max-amplicon-size 460 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
  --R1-size 250 --R2-size 250 --expected-amplicon-size 420 \
 --input-archive $frogs_dir/test/data/test_dataset.tar.gz \
 --output-dereplicated $frogs_dir/test/results/01-prepro.fasta \
 --output-count $frogs_dir/test/results/01-prepro.tsv \
 --summary $frogs_dir/test/results/01-prepro.html \
 --log-file $frogs_dir/test/results/01-prepro.log \
 --nb-cpus $nb_cpu 
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess" >&2
	exit 1;
fi


echo "Step clustering `date`"

$frogs_dir/tools/clustering.py \
 --distance 3 \
 --denoising \
 --input-fasta $frogs_dir/test/results/01-prepro.fasta \
 --input-count $frogs_dir/test/results/01-prepro.tsv \
 --output-biom $frogs_dir/test/results/02-clustering.biom \
 --output-fasta $frogs_dir/test/results/02-clustering.fasta \
 --output-compo $frogs_dir/test/results/02-clustering_compo.tsv \
 --log-file $frogs_dir/test/results/02-clustering.log \
 --nb-cpus $nb_cpu

if [ $? -ne 0 ]
then
	echo "Error in clustering" >&2
	exit 1;
fi


echo "Step remove_chimera `date`"

$frogs_dir/tools/remove_chimera.py \
 --input-fasta $frogs_dir/test/results/02-clustering.fasta \
 --input-biom $frogs_dir/test/results/02-clustering.biom \
 --non-chimera $frogs_dir/test/results/03-chimera.fasta \
 --out-abundance $frogs_dir/test/results/03-chimera.biom \
 --summary $frogs_dir/test/results/03-chimera.html \
 --log-file $frogs_dir/test/results/03-chimera.log \
 --nb-cpus $nb_cpu
 
if [ $? -ne 0 ]
then
	echo "Error in remove_chimera" >&2
	exit 1;
fi


echo "Step filters `date`"

$frogs_dir/tools/filters.py \
 --min-abundance 0.00005 \
 --input-biom $frogs_dir/test/results/03-chimera.biom \
 --input-fasta $frogs_dir/test/results/03-chimera.fasta \
 --output-fasta $frogs_dir/test/results/04-filters.fasta \
 --output-biom $frogs_dir/test/results/04-filters.biom \
 --excluded $frogs_dir/test/results/04-filters.excluded \
 --summary $frogs_dir/test/results/04-filters.html \
 --log-file $frogs_dir/test/results/04-filters.log \

if [ $? -ne 0 ]
then
	echo "Error in filters" >&2
	exit 1;
fi


echo "Step affiliation_OTU `date`"

$frogs_dir/tools/affiliation_OTU.py \
 --reference $frogs_dir/test/data/db.fasta \
 --input-fasta $frogs_dir/test/results/04-filters.fasta \
 --input-biom $frogs_dir/test/results/04-filters.biom \
 --output-biom $frogs_dir/test/results/04-affiliation.biom \
 --summary $frogs_dir/test/results/04-affiliation.html \
 --log-file $frogs_dir/test/results/04-affiliation.log \
 --nb-cpus $nb_cpu --java-mem $java_mem

if [ $? -ne 0 ]
then
	echo "Error in affiliation_OTU" >&2
	exit 1;
fi


echo "Step clusters_stat `date`"

$frogs_dir/tools/clusters_stat.py \
 --input-biom $frogs_dir/test/results/04-affiliation.biom \
 --output-file $frogs_dir/test/results/05-clustersStat.html \
 --log-file $frogs_dir/test/results/05-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stat" >&2
	exit 1;
fi


echo "Step affiliations_stat `date`"

$frogs_dir/tools/affiliations_stat.py \
 --input-biom $frogs_dir/test/results/04-affiliation.biom \
 --output-file $frogs_dir/test/results/06-affiliationsStat.html \
 --log-file $frogs_dir/test/results/06-affiliationsStat.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species

if [ $? -ne 0 ]
then
	echo "Error in affiliations_stat" >&2
	exit 1;
fi


echo "Step biom_to_tsv `date`"

$frogs_dir/tools/biom_to_tsv.py \
 --input-biom $frogs_dir/test/results/04-affiliation.biom \
 --input-fasta $frogs_dir/test/results/04-filters.fasta \
 --output-tsv $frogs_dir/test/results/07-biom2tsv.tsv \
 --output-multi-affi $frogs_dir/test/results/07-biom2tsv.multi \
 --log-file $frogs_dir/test/results/07-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi


echo "Completed with success"
