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
export PATH=$frogs_dir/libexec:$frogs_dir/app:$PATH
export PYTHONPATH=$frogs_dir/lib:$PYTHONPATH


# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi


echo "Step preprocess : Flash `date`"

preprocess.py illumina \
 --min-amplicon-size 44 --max-amplicon-size 490 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
 --R1-size 267 --R2-size 266 --expected-amplicon-size 420 --merge-software flash\
 --nb-cpus $nb_cpu --mismatch-rate 0.15 --keep-unmerged \
 --input-archive $frogs_dir/test/data/test_dataset.tar.gz \
 --output-dereplicated $out_dir/01-prepro-flash.fasta \
 --output-count $out_dir/01-prepro-flash.tsv \
 --summary $out_dir/01-prepro-flash.html \
 --log-file $out_dir/01-prepro-flash.log 
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess : Flash" >&2
	exit 1;
fi

echo "Step preprocess : Vsearch `date`"

preprocess.py illumina \
 --min-amplicon-size 44 --max-amplicon-size 490 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
 --R1-size 267 --R2-size 266 --expected-amplicon-size 420 --merge-software flash\
 --nb-cpus $nb_cpu --mismatch-rate 0.15 --keep-unmerged \
 --input-archive $frogs_dir/test/data/test_dataset.tar.gz \
 --output-dereplicated $out_dir/01-prepro-vsearch.fasta \
 --output-count $out_dir/01-prepro-vsearch.tsv \
 --summary $out_dir/01-prepro-vsearch.html \
 --log-file $out_dir/01-prepro-vsearch.log 
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess : Vsearch" >&2
	exit 1;
fi

echo "Step clustering `date`"

clustering.py \
 --distance 3 \
 --denoising \
 --input-fasta $out_dir/01-prepro-vsearch.fasta \
 --input-count $out_dir/01-prepro-vsearch.tsv \
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

remove_chimera.py \
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

filters.py \
 --min-abundance 0.00005 \
 --min-sample-presence 3 \
 --nb-cpus $nb_cpu \
 --input-biom $out_dir/03-chimera.biom \
 --input-fasta $out_dir/03-chimera.fasta \
 --output-fasta $out_dir/04-filters.fasta \
 --output-biom $out_dir/04-filters.biom \
 --excluded $out_dir/04-filters.excluded \
 --summary $out_dir/04-filters.html \
 --log-file $out_dir/04-filters.log 

if [ $? -ne 0 ]
then
	echo "Error in filters" >&2
	exit 1;
fi

echo "Step ITSx `date`"

itsx.py \
 --input-fasta $out_dir/04-filters.fasta \
 --input-biom $out_dir/04-filters.biom \
 --region ITS1 --nb-cpus $nb_cpu \
 --out-abundance $out_dir/05-itsx.biom \
 --summary $out_dir/05-itsx.html \
 --log-file $out_dir/05-itsx.log \
 --out-fasta $out_dir/05-itsx.fasta \
 --out-removed $out_dir/05-itsx-excluded.tsv

if [ $? -ne 0 ]
then
	echo "Error in ITSx" >&2
	exit 1;
fi

echo "Step affiliation_OTU `date`"

affiliation_OTU.py \
 --reference $frogs_dir/test/data/ITS1.rdp.fasta \
 --input-fasta $out_dir/04-filters.fasta \
 --input-biom $out_dir/04-filters.biom \
 --output-biom $out_dir/06-affiliation.biom \
 --summary $out_dir/06-affiliation.html \
 --log-file $out_dir/06-affiliation.log \
 --nb-cpus $nb_cpu --java-mem $java_mem \
 --rdp

if [ $? -ne 0 ]
then
	echo "Error in affiliation_OTU" >&2
	exit 1;
fi

echo "Step affiliation_postprocess `date`"

affiliation_postprocess.py \
 --input-biom $out_dir/06-affiliation.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --reference $frogs_dir/test/data/Unite_extract_ITS1.fasta \
 --output-biom $out_dir/07-affiliation_postprocessed.biom \
 --output-compo $out_dir/07-affiliation_postprocessed.compo.tsv \
 --output-fasta $out_dir/07-affiliation_postprocessed.fasta \
 --log-file $out_dir/07-affiliation_postprocessed.log

if [ $? -ne 0 ]
then
	echo "Error in affiliation_postprocess" >&2
	exit 1;
fi


echo "Step clusters_stat `date`"

clusters_stat.py \
 --input-biom $out_dir/06-affiliation.biom \
 --output-file $out_dir/08-clustersStat.html \
 --log-file $out_dir/08-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stat" >&2
	exit 1;
fi


echo "Step affiliations_stat `date`"

affiliations_stat.py \
 --input-biom $out_dir/06-affiliation.biom \
 --output-file $out_dir/09-affiliationsStat.html \
 --log-file $out_dir/09-affiliationsStat.log \
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

biom_to_tsv.py \
 --input-biom $out_dir/06-affiliation.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --output-tsv $out_dir/10-biom2tsv.tsv \
 --output-multi-affi $out_dir/10-biom2tsv-affiliation_multihit.tsv \
 --log-file $out_dir/10-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi

echo "Step biom_to_stdBiom `date`"


biom_to_stdBiom.py \
 --input-biom $out_dir/06-affiliation.biom \
 --output-biom $out_dir/11-affiliation_std.biom \
 --output-metadata $out_dir/11-affiliation_multihit.tsv \
 --log-file $out_dir/11-biom2stdbiom.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_stdBiom" >&2
	exit 1;
fi

echo "Step tsv_to_biom `date`"


tsv_to_biom.py \
 --input-tsv $out_dir/10-biom2tsv.tsv \
 --input-multi-affi $out_dir/10-biom2tsv-affiliation_multihit.tsv \
 --output-biom $out_dir/12-tsv2biom.biom \
 --output-fasta $out_dir/12-tsv2biom.fasta \
 --log-file $out_dir/12-tsv2biom.log 

if [ $? -ne 0 ]
then
	echo "Error in tsv_to_biom" >&2
	exit 1;
fi

echo "Step tree : pynast `date`"

tree.py \
 --nb-cpus $nb_cpu  \
 --input-otu $frogs_dir/test/data/to_test_pynast.fasta \
 --biomfile $frogs_dir/test/data/to_test_pynast.biom \
 --template-pynast $frogs_dir/test/data/otus_pynast.fasta \
 --out-tree $out_dir/13-tree-pynast.nwk \
 --html $out_dir/13-tree-pynast.html \
 --log-file $out_dir/13-tree-pynast.log
 
if [ $? -ne 0 ]
then
	echo "Error in tree : pynast" >&2
	exit 1;
fi

echo "Step tree : mafft `date`"

tree.py \
 --nb-cpus $nb_cpu \
 --input-otu $out_dir/04-filters.fasta \
 --biomfile $out_dir/06-affiliation.biom \
 --out-tree $out_dir/13-tree-mafft.nwk \
 --html $out_dir/13-tree-mafft.html \
 --log-file $out_dir/13-tree-mafft.log

if [ $? -ne 0 ]
then
	echo "Error in tree : mafft" >&2
	exit 1;
fi

echo "Step r_import_data `date`"

r_import_data.py  \
 --biomfile $frogs_dir/test/data/chaillou.biom \
 --samplefile $frogs_dir/test/data/sample_metadata.tsv \
 --treefile $frogs_dir/test/data/tree.nwk \
 --rdata $out_dir/14-phylo_import.Rdata --html $out_dir/14-phylo_import.html --log-file $out_dir/14-phylo_import.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_import_data " >&2
	exit 1;
fi

echo "Step r_composition `date`"

r_composition.py  \
 --varExp EnvType --taxaRank1 Kingdom --taxaSet1 Bacteria --taxaRank2 Phylum --numberOfTaxa 9 \
 --rdata $out_dir/14-phylo_import.Rdata \
 --html $out_dir/15-phylo_composition.html --log-file $out_dir/15-phylo_composition.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_composition " >&2
	exit 1;
fi

echo "Step r_alpha_diversity `date`"

r_alpha_diversity.py  \
 --varExp EnvType \
 --rdata $out_dir/14-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --alpha-out $out_dir/16-phylo_alpha_div.tsv --html $out_dir/16-phylo_alpha_div.html --log-file $out_dir/16-phylo_alpha_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_alpha_diversity " >&2
	exit 1;
fi

echo "Step r_beta_diversity `date`"

r_beta_diversity.py  \
 --varExp EnvType --distance-methods cc,unifrac \
 --rdata $out_dir/14-phylo_import.Rdata \
 --matrix-outdir $out_dir --html $out_dir/17-phylo_beta_div.html --log-file $out_dir/17-phylo_beta_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_beta_diversity " >&2
	exit 1;
fi

echo "Step r_structure `date`"

r_structure.py  \
 --varExp EnvType --ordination-method MDS \
 --rdata $out_dir/14-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/18-phylo_structure.html --log-file $out_dir/18-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_structure " >&2
	exit 1;
fi

echo "Step r_clustering `date`"

r_clustering.py  \
 --varExp EnvType \
 --rdata $out_dir/14-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/19-phylo_clutering.html --log-file $out_dir/16-phylo_clustering.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_clustering " >&2
	exit 1;
fi

echo "Step r_manova `date`"

r_manova.py  \
 --varExp EnvType \
 --rdata $out_dir/14-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/20-phylo_manova.html --log-file $out_dir/20-phylo_manova.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_manova " >&2
	exit 1;
fi

echo "Completed with success"
