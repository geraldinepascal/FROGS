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
 --R1-size 267 --R2-size 266 --expected-amplicon-size 420 --merge-software flash \
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
 --R1-size 267 --R2-size 266 --merge-software vsearch \
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


echo "Step otu filters `date`"

otu-filters.py \
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


echo "Step affiliation_filter: masking mode `date`"

affiliation_filters.py \
--input-biom $out_dir/06-affiliation.biom \
--input-fasta $out_dir/04-filters.fasta \
--output-biom $out_dir/07-affiliation_masked.biom \
--summary $out_dir/07-affiliation_masked.html \
--impacted $out_dir/07-impacted_OTU_masked.tsv \
--log-file $out_dir/07-affiliation_filter_maskMode.log \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 150 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 1e-150 \
--taxon-ignored "g__Sarcodon" "s__Trichoderma" \
--mask \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species  # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_filter" >&2
	exit 1;
fi


echo "Step affiliation_filter: deleted mode `date`"

affiliation_filters.py \
--input-biom $out_dir/06-affiliation.biom \
--input-fasta $out_dir/04-filters.fasta \
--output-biom $out_dir/07-affiliation_deleted.biom \
--output-fasta $out_dir/07-affiliation_deleted.fasta \
--summary $out_dir/07-affiliation_deleted.html \
--impacted $out_dir/07-impacted_OTU_deleted.tsv \
--log-file $out_dir/07-affiliation_filter_delMode.log \
--min-rdp-bootstrap Species:0.8 \
--min-blast-length 150 \
--min-blast-identity 1.0 \
--min-blast-coverage 1.0 \
--max-blast-evalue 1e-150 \
--taxon-ignored "g__Sarcodon" "s__Trichoderma" \
--delete \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species  # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliation_filter" >&2
	exit 1;
fi


echo "Step affiliation_postprocess `date`"

affiliation_postprocess.py \
 --input-biom $out_dir/06-affiliation.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --reference $frogs_dir/test/data/Unite_extract_ITS1.fasta \
 --output-biom $out_dir/08-affiliation_postprocessed.biom \
 --output-compo $out_dir/08-affiliation_postprocessed.compo.tsv \
 --output-fasta $out_dir/08-affiliation_postprocessed.fasta \
 --log-file $out_dir/08-affiliation_postprocessed.log

if [ $? -ne 0 ]
then
	echo "Error in affiliation_postprocess" >&2
	exit 1;
fi


echo "Step clusters_stat `date`"

clusters_stat.py \
 --input-biom $out_dir/07-affiliation_masked.biom \
 --output-file $out_dir/09-clustersStat.html \
 --log-file $out_dir/09-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stat" >&2
	exit 1;
fi

					
echo "Step affiliations_stat `date`"

affiliations_stat.py \
 --input-biom $out_dir/07-affiliation_masked.biom \
 --output-file $out_dir/10-affiliationsStat.html \
 --log-file $out_dir/10-affiliationsStat.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species \
 --taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

if [ $? -ne 0 ]
then
	echo "Error in affiliations_stat" >&2
	exit 1;
fi


echo "Step biom_to_tsv `date`"

biom_to_tsv.py \
 --input-biom $out_dir/07-affiliation_masked.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --output-tsv $out_dir/11-biom2tsv.tsv \
 --output-multi-affi $out_dir/11-biom2tsv-affiliation_multihit.tsv \
 --log-file $out_dir/11-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi

echo "Step biom_to_stdBiom `date`"


biom_to_stdBiom.py \
 --input-biom $out_dir/07-affiliation_masked.biom \
 --output-biom $out_dir/12-affiliation_std.biom \
 --output-metadata $out_dir/12-affiliation_multihit.tsv \
 --log-file $out_dir/12-biom2stdbiom.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_stdBiom" >&2
	exit 1;
fi

echo "Step tsv_to_biom `date`"


tsv_to_biom.py \
 --input-tsv $out_dir/11-biom2tsv.tsv \
 --input-multi-affi $out_dir/11-biom2tsv-affiliation_multihit.tsv \
 --output-biom $out_dir/13-tsv2biom.biom \
 --output-fasta $out_dir/13-tsv2biom.fasta \
 --log-file $out_dir/13-tsv2biom.log 

if [ $? -ne 0 ]
then
	echo "Error in tsv_to_biom" >&2
	exit 1;
fi

echo "Step tree : mafft `date`"

tree.py \
 --nb-cpus $nb_cpu \
 --input-sequences $out_dir/04-filters.fasta \
 --biom-file $out_dir/07-affiliation_masked.biom \
 --out-tree $out_dir/14-tree-mafft.nwk \
 --html $out_dir/14-tree-mafft.html \
 --log-file $out_dir/14-tree-mafft.log

if [ $? -ne 0 ]
then
	echo "Error in tree : mafft" >&2
	exit 1;
fi

echo "Step phyloseq_import_data `date`"

phyloseq_import_data.py  \
 --biomfile $frogs_dir/test/data/chaillou.biom \
 --samplefile $frogs_dir/test/data/sample_metadata.tsv \
 --treefile $frogs_dir/test/data/tree.nwk \
 --rdata $out_dir/15-phylo_import.Rdata \
 --html $out_dir/15-phylo_import.html \
 --log-file $out_dir/15-phylo_import.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_import_data " >&2
	exit 1;
fi

echo "Step phyloseq_composition `date`"

phyloseq_composition.py  \
 --varExp EnvType --taxaRank1 Kingdom --taxaSet1 Bacteria --taxaRank2 Phylum --numberOfTaxa 9 \
 --rdata $out_dir/15-phylo_import.Rdata \
 --html $out_dir/16-phylo_composition.html \
 --log-file $out_dir/16-phylo_composition.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_composition " >&2
	exit 1;
fi

echo "Step phyloseq_alpha_diversity `date`"

phyloseq_alpha_diversity.py  \
 --varExp EnvType \
 --rdata $out_dir/15-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --alpha-out $out_dir/17-phylo_alpha_div.tsv \
 --html $out_dir/17-phylo_alpha_div.html \
 --log-file $out_dir/17-phylo_alpha_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_alpha_diversity " >&2
	exit 1;
fi

echo "Step phyloseq_beta_diversity `date`"

phyloseq_beta_diversity.py  \
 --varExp EnvType --distance-methods cc,unifrac \
 --rdata $out_dir/15-phylo_import.Rdata \
 --matrix-outdir $out_dir \
 --html $out_dir/18-phylo_beta_div.html \
 --log-file $out_dir/18-phylo_beta_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_beta_diversity " >&2
	exit 1;
fi

echo "Step phyloseq_structure `date`"

phyloseq_structure.py  \
 --varExp EnvType --ordination-method MDS \
 --rdata $out_dir/15-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/19-phylo_structure.html \
 --log-file $out_dir/19-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_structure " >&2
	exit 1;
fi

echo "Step phyloseq_clustering `date`"

phyloseq_clustering.py  \
 --varExp EnvType \
 --rdata $out_dir/15-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/20-phylo_clutering.html \
 --log-file $out_dir/20-phylo_clustering.log

 
if [ $? -ne 0 ]
then
	echo "Error in phyloseq_clustering " >&2
	exit 1;
fi

echo "Step phyloseq_manova `date`"

phyloseq_manova.py  \
 --varExp EnvType \
 --rdata $out_dir/15-phylo_import.Rdata --distance-matrix $out_dir/unifrac.tsv \
 --html $out_dir/21-phylo_manova.html \
 --log-file $out_dir/21-phylo_manova.log

if [ $? -ne 0 ]
then
	echo "Error in phyloseq_manova " >&2
	exit 1;
fi


echo "Step deseq2_preprocess `date`"

deseq2_preprocess.py \
 --data $out_dir/15-phylo_import.Rdata \
 --log-file $out_dir/22-deseq2_preprocess.log \
 --out-Rdata $out_dir/22-deseq2_preprocess.Rdata \
 --var EnvType

if [ $? -ne 0 ]
then
	echo "Error in deseq2_preprocess " >&2
	exit 1;
fi


echo "Step deseq2_visualization `date`"

deseq2_visualization.py \
 --phyloseqData $out_dir/15-phylo_import.Rdata \
 --dds $out_dir/22-deseq2_preprocess.Rdata \
 --log-file $out_dir/23-deseq2_visualization.log \
 --html $out_dir/22-deseq2_visualization.html \
 --var EnvType --mod1 BoeufHache --mod2 SaumonFume
                            

if [ $? -ne 0 ]
then
	echo "Error in deseq2_visualization " >&2
	exit 1;
fi

echo "Step normalisation `date`"
./normalisation.py \
 -n 100 \
 -i $out_dir/07-affiliation_masked.biom \
 -f $out_dir/07-affiliation_deleted.fasta \
 -b $out_dir/test_aff_N1000.biom \
 -o $out_dir/24-test_aff_N1000.fasta \
 -s $out_dir/24-normalisation.html \
 -l $out_dir/24-normalisation.log
 
if [ $? -ne 0 ]
then
	echo "Error in normalisation" >&2
	exit 1;
fi

echo "Completed with success"

