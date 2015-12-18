#!/bin/bash
frogs_dir=$1
nb_cpu=$2
java_mem=$3
out_dir=$4

# Check parameters
if [ "$#" -ne 4 ]; then
    echo "ERROR: Illegal number of parameters." ;
    echo 'Command usage: test.sh <FROGS_FOLDER> <NB_CPU> <JAVA_MEM> <OUT_FOLDER>' ;
    exit 1 ;
fi

# Set ENV
export PATH=$frogs_dir/bin:$PATH
export PYTHONPATH=$frogs_dir/bin:$PYTHONPATH


# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi


echo "Step preprocess `date`"

$frogs_dir/tools/preprocess/preprocess.py illumina \
 --min-amplicon-size 380 --max-amplicon-size 460 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
  --R1-size 250 --R2-size 250 --expected-amplicon-size 420 \
 --input-archive $frogs_dir/test/data/test_dataset.tar.gz \
 --output-dereplicated $out_dir/01-prepro.fasta \
 --output-count $out_dir/01-prepro.tsv \
 --summary $out_dir/01-prepro.html \
 --log-file $out_dir/01-prepro.log \
 --nb-cpus $nb_cpu 
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess" >&2
	exit 1;
fi


echo "Step clustering `date`"

$frogs_dir/tools/clustering/clustering.py \
 --distance 3 \
 --denoising \
 --input-fasta $out_dir/01-prepro.fasta \
 --input-count $out_dir/01-prepro.tsv \
 --output-biom $out_dir/02-clustering.biom \
 --output-fasta $out_dir/02-clustering.fasta \
 --output-compo $out_dir/02-clustering_compo.tsv \
 --log-file $out_dir/02-clustering.log \
 --nb-cpus $nb_cpu

if [ $? -ne 0 ]
then
	echo "Error in clustering" >&2
	exit 1;
fi


echo "Step remove_chimera `date`"

$frogs_dir/tools/remove_chimera/remove_chimera.py \
 --input-fasta $out_dir/02-clustering.fasta \
 --input-biom $out_dir/02-clustering.biom \
 --non-chimera $out_dir/03-chimera.fasta \
 --out-abundance $out_dir/03-chimera.biom \
 --summary $out_dir/03-chimera.html \
 --log-file $out_dir/03-chimera.log \
 --nb-cpus $nb_cpu
 
if [ $? -ne 0 ]
then
	echo "Error in remove_chimera" >&2
	exit 1;
fi


echo "Step filters `date`"

$frogs_dir/tools/filters/filters.py \
 --min-abundance 0.00005 \
 --input-biom $out_dir/03-chimera.biom \
 --input-fasta $out_dir/03-chimera.fasta \
 --output-fasta $out_dir/04-filters.fasta \
 --output-biom $out_dir/04-filters.biom \
 --excluded $out_dir/04-filters.excluded \
 --summary $out_dir/04-filters.html \
 --log-file $out_dir/04-filters.log \

if [ $? -ne 0 ]
then
	echo "Error in filters" >&2
	exit 1;
fi


echo "Step affiliation_OTU `date`"

$frogs_dir/tools/affiliation_OTU/affiliation_OTU.py \
 --reference $frogs_dir/test/data/db.fasta \
 --input-fasta $out_dir/04-filters.fasta \
 --input-biom $out_dir/04-filters.biom \
 --output-biom $out_dir/04-affiliation.biom \
 --summary $out_dir/04-affiliation.html \
 --log-file $out_dir/04-affiliation.log \
 --nb-cpus $nb_cpu --java-mem $java_mem

if [ $? -ne 0 ]
then
	echo "Error in affiliation_OTU" >&2
	exit 1;
fi


echo "Step clusters_stat `date`"

$frogs_dir/tools/clusters_stat/clusters_stat.py \
 --input-biom $out_dir/04-affiliation.biom \
 --output-file $out_dir/05-clustersStat.html \
 --log-file $out_dir/05-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stat" >&2
	exit 1;
fi


echo "Step affiliations_stat `date`"

$frogs_dir/tools/affiliations_stat/affiliations_stat.py \
 --input-biom $out_dir/04-affiliation.biom \
 --output-file $out_dir/06-affiliationsStat.html \
 --log-file $out_dir/06-affiliationsStat.log \
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

$frogs_dir/tools/biom_to_tsv/biom_to_tsv.py \
 --input-biom $out_dir/04-affiliation.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --output-tsv $out_dir/07-biom2tsv.tsv \
 --output-multi-affi $out_dir/07-biom2tsv.multi \
 --log-file $out_dir/07-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi


echo "Completed with success"
