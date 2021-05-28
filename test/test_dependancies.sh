#!/bin/bash

nb_cpu=2
java_mem=1
out_dir=res_3.2.3_ter
expected_dir=res_3.2.3
run_programs=true	## if true lance les python sinon, fait uniquement les comparatifs de résultats

## Set ENV
## export PATH=../app:$PATH
## or
## conda activate __frogs@3.2.0
##		# après installation de l'env, check de la présence des librairies perl

## Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi

diff_line() {
	sort $1 > /tmp/1_`basename $1`
	sort $2 > /tmp/2_`basename $2`
	n_diff=`sdiff -s /tmp/1_$(basename $1) /tmp/2_$(basename $2) | wc -l`
	rm /tmp/1_`basename $1` /tmp/2_`basename $2`
	if [ "$n_diff" -gt $3  ]; then true; else false; fi
}

diff_size() {
  	diff=`ls -l $1 $2 | awk '{if(a==""){a=$5}else{b=$5}}END{if(a<=b){print b-a}else{print a-b}}' `
	if [ "$diff" -gt $3  ];	then true; else false; fi
}

echo "Step demultiplexe `date`"
demultiplex.py \
  --input-R1 data/demultiplex.fastq.gz --input-barcode data/demultiplex.barcode.txt \
  --mismatches 1 --end both \
  --output-demultiplexed $out_dir/demultiplexed.tar.gz --output-excluded $out_dir/undemultiplexed.tar.gz --log-file $out_dir/demultiplex.log --summary $out_dir/demultiplex_summary.txt 

if diff_line $out_dir/demultiplex_summary.txt $expected_dir/demultiplex_summary.txt  0
then
	echo "difference in clustering: 02-clustering_fastidious.biom " >&2
fi

echo "Step preprocess : Flash `date`"

if $run_programs
then
	preprocess.py illumina \
	 --min-amplicon-size 44 --max-amplicon-size 490 \
	 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
	 --R1-size 267 --R2-size 266 --expected-amplicon-size 420 --merge-software flash \
	 --nb-cpus $nb_cpu --mismatch-rate 0.15 --keep-unmerged \
	 --input-archive data/test_dataset.tar.gz \
	 --output-dereplicated $out_dir/01-prepro-flash.fasta \
	 --output-count $out_dir/01-prepro-flash.tsv \
	 --summary $out_dir/01-prepro-flash.html \
	 --log-file $out_dir/01-prepro-flash.log
	 
	if [ $? -ne 0 ]
	then
		echo "Error in preprocess : Flash" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/01-prepro-flash.fasta $expected_dir/01-prepro-flash.fasta 0
then
	echo "difference in preprocess, Flash: 01-prepro-flash.fasta" >&2
fi

if diff_line $out_dir/01-prepro-flash.tsv $expected_dir/01-prepro-flash.tsv 0
then
	echo "difference in preprocess, Flash: 01-prepro-flash.tsv " >&2
fi

echo "Step preprocess : Vsearch `date`"

if $run_programs
then
	preprocess.py illumina \
	 --min-amplicon-size 44 --max-amplicon-size 490 \
	 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
	 --R1-size 267 --R2-size 266 --merge-software vsearch \
	 --nb-cpus $nb_cpu --mismatch-rate 0.15 --keep-unmerged \
	 --input-archive data/test_dataset.tar.gz \
	 --output-dereplicated $out_dir/01-prepro-vsearch.fasta \
	 --output-count $out_dir/01-prepro-vsearch.tsv \
	 --summary $out_dir/01-prepro-vsearch.html \
	 --log-file $out_dir/01-prepro-vsearch.log 
	 
	if [ $? -ne 0 ]
	then
		echo "Error in preprocess : Vsearch" >&2
		exit 1;
	fi
fi

 
if diff_line $out_dir/01-prepro-vsearch.fasta $expected_dir/01-prepro-vsearch.fasta 0
then
	echo "difference in preprocess, Flash: 01-prepro-vsearch.fasta" >&2
fi

if diff_line $out_dir/01-prepro-vsearch.tsv $expected_dir/01-prepro-vsearch.tsv 0
then
	echo "difference in preprocess, Flash: 01-prepro-vsearch.tsv " >&2
fi

if diff_line $out_dir/01-prepro-vsearch.html $expected_dir/01-prepro-vsearch.html 0
then
	echo "difference in preprocess, Flash: 01-prepro-vsearch.html " >&2
fi

echo "Step clustering fastidious `date`"

if $run_programs
then
	clustering.py \
	 --distance 1 \
	 --fastidious \
	 --input-fasta $expected_dir/01-prepro-vsearch.fasta \
	 --input-count $expected_dir/01-prepro-vsearch.tsv \
	 --output-biom $out_dir/02-clustering_fastidious.biom \
	 --output-fasta $out_dir/02-clustering_fastidious.fasta \
	 --output-compo $out_dir/02-clustering_fastidious_compo.tsv \
	 --log-file $out_dir/02-clustering_fastidious.log \
	 --nb-cpus $nb_cpu

	if [ $? -ne 0 ]
	then
		echo "Error in clustering fastidious" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/02-clustering_fastidious.fasta $expected_dir/02-clustering_fastidious.fasta 0
then
	echo "difference in clustering: 02-clustering_fastidious.fasta " >&2
fi

if diff_size $out_dir/02-clustering_fastidious.biom $expected_dir/02-clustering_fastidious.biom 1
then
	echo "difference in clustering: 02-clustering_fastidious.biom " >&2
fi

if diff_line $out_dir/02-clustering_fastidious_compo.tsv $expected_dir/02-clustering_fastidious_compo.tsv 0
then
	echo "difference in clustering: 02-clustering_fastidious_compo.tsv " >&2
fi

echo "Step remove_chimera `date`"

if $run_programs
then
	remove_chimera.py \
	 --input-fasta $expected_dir/02-clustering_fastidious.fasta \
	 --input-biom $expected_dir/02-clustering_fastidious.biom \
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
fi

if diff_line $out_dir/03-chimera.fasta $expected_dir/03-chimera.fasta 0
then
	echo "difference in remove_chimera: 03-chimera.fasta " >&2
fi

if diff_size $out_dir/03-chimera.biom $expected_dir/03-chimera.biom 0
then
	echo "difference in remove_chimera: 03-chimera.biom " >&2
fi

if diff_line $out_dir/03-chimera.html $expected_dir/03-chimera.html 0
then
	echo "difference in remove_chimera: 03-chimera.html " >&2
fi

echo "Step otu filters `date`"

if $run_programs
then
	otu_filters.py \
	 --min-abundance 0.00005 \
	 --min-sample-presence 3 \
	 --contaminant data/phi.fa \
	 --nb-cpus $nb_cpu \
	 --input-biom $expected_dir/03-chimera.biom \
	 --input-fasta $expected_dir/03-chimera.fasta \
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
fi

if diff_line $out_dir/04-filters.fasta $expected_dir/04-filters.fasta 0
then
	echo "difference in otu_filters : 04-filters.fasta " >&2
fi

if diff_size $out_dir/04-filters.biom $expected_dir/04-filters.biom 0
then
	echo "difference in otu_filters : 04-filters.biom" >&2
fi

if diff_line $out_dir/04-filters.excluded $expected_dir/04-filters.excluded 0
then
	echo "difference in otu_filters : 04-filters.excluded " >&2
fi

if diff_line $out_dir/04-filters.html $expected_dir/04-filters.html 0
then
	echo "difference in otu_filters : 04-filters.html " >&2
fi

echo "Step ITSx `date`"
if $run_programs
then
	itsx.py \
	 --input-fasta $expected_dir/04-filters.fasta \
	 --input-biom $expected_dir/04-filters.biom \
	 --region ITS1 --check-its-only --nb-cpus $nb_cpu \
	 --out-abundance $out_dir/05-itsx.biom \
	 --summary $out_dir/05-itsx.html \
	 --log-file $out_dir/05-itsx.log \
	 --out-fasta $out_dir/05-itsx.fasta \
	 --out-removed $out_dir/05-itsx-excluded.fasta

	if [ $? -ne 0 ]
	then
		echo "Error in ITSx" >&2
		exit 1;
	fi
fi

if diff_size $out_dir/05-itsx.biom $expected_dir/05-itsx.biom 0
then
	echo "difference in ITSx : 05-itsx.biom " >&2
fi

if diff_line $out_dir/05-itsx.fasta $expected_dir/05-itsx.fasta 0
then
	echo "difference in ITSx : 05-itsx.fasta " >&2
fi

if diff_line $out_dir/05-itsx-excluded.fasta $expected_dir/05-itsx-excluded.fasta 0
then
	echo "difference in ITSx : 05-itsx-excluded.fasta " >&2
fi

if diff_line $out_dir/05-itsx.html $expected_dir/05-itsx.html 0
then
	echo "difference in ITSx : 05-itsx.html " >&2
fi

echo "Step affiliation_OTU `date`"

if $run_programs
then
	affiliation_OTU.py \
	 --reference data/ITS1.rdp.fasta \
	 --input-fasta $expected_dir/04-filters.fasta \
	 --input-biom $expected_dir/04-filters.biom \
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
fi

if diff_size $out_dir/06-affiliation.biom $expected_dir/06-affiliation.biom 0
then
	echo "difference in affiliation_OTU :06-affiliation.biom " >&2
fi

if diff_line $out_dir/06-affiliation.html $expected_dir/06-affiliation.html 0
then
	echo "difference in affiliation_OTU :06-affiliation.html " >&2
fi

echo "Step affiliation_filter: masking mode `date`"

if $run_programs
then
	affiliation_filters.py \
	--input-biom $expected_dir/06-affiliation.biom \
	--input-fasta $expected_dir/04-filters.fasta \
	--output-biom $out_dir/07-affiliation_masked.biom \
	--summary $out_dir/07-affiliation_masked.html \
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
	--taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

	if [ $? -ne 0 ]
	then
		echo "Error in affiliation_filter" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/07-impacted_OTU_masked.tsv $expected_dir/07-impacted_OTU_masked.tsv 0
then
	echo "difference in affiliation_filters, mask mode :07-impacted_OTU_masked.tsv" >&2
fi

if diff_line $out_dir/07-affiliation_masked.html $expected_dir/07-affiliation_masked.html 0
then
	echo "difference in affiliation_filters, mask mode :07-affiliation_masked.html" >&2
fi

if diff_line $out_dir/07-impacted_OTU_masked_multihit.tsv $expected_dir/07-impacted_OTU_masked_multihit.tsv 0
then
	echo "difference in affiliation_filters, mask mode :07-impacted_OTU_masked_multihit.tsv" >&2
fi

if diff_size $out_dir/07-affiliation_masked.biom  $expected_dir/07-affiliation_masked.biom 0
then
	echo "difference in affiliation_filters, mask mode :07-affiliation_masked.biom" >&2
fi

echo "Step affiliation_filter: deleted mode `date`"

if $run_programs
then
	affiliation_filters.py \
	--input-biom $expected_dir/06-affiliation.biom \
	--input-fasta $expected_dir/04-filters.fasta \
	--output-biom $out_dir/07-affiliation_deleted.biom \
	--output-fasta $out_dir/07-affiliation_deleted.fasta \
	--summary $out_dir/07-affiliation_deleted.html \
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
	--taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

	if [ $? -ne 0 ]
	then
		echo "Error in affiliation_filter" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/07-impacted_OTU_deleted.tsv $expected_dir/07-impacted_OTU_deleted.tsv 0
then
	echo "difference in affiliation_filters, delete mode :07-impacted_OTU_deleted.tsv" >&2
fi

if diff_line $out_dir/07-affiliation_deleted.html $expected_dir/07-affiliation_deleted.html 0
then
	echo "difference in affiliation_filters, delete mode :07-affiliation_deleted.html" >&2
fi

if diff_line $out_dir/07-impacted_OTU_deleted_multihit.tsv $expected_dir/07-impacted_OTU_deleted_multihit.tsv 0
then
	echo "difference in affiliation_filters, delete mode :07-impacted_OTU_deleted_multihit.tsv" >&2
fi

if diff_size $out_dir/07-affiliation_deleted.biom  $expected_dir/07-affiliation_deleted.biom 0
then
	echo "difference in affiliation_filters, delete mode :07-affiliation_deleted.biom" >&2
fi

if diff_line $out_dir/07-affiliation_deleted.fasta $expected_dir/07-affiliation_deleted.fasta 0
then
	echo "difference in affiliation_filters, delete mode :$out_dir/07-affiliation_deleted.fasta" >&2
fi

echo "Step affiliation_postprocess `date`"

if $run_programs
then
	affiliation_postprocess.py \
	 --input-biom $expected_dir/06-affiliation.biom \
	 --input-fasta $expected_dir/04-filters.fasta \
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
fi

if diff_line $out_dir/08-affiliation_postprocessed.compo.tsv $expected_dir/08-affiliation_postprocessed.compo.tsv 0
then
	echo "difference in affiliation_postprocess :08-affiliation_postprocessed.compo.tsv" >&2
fi	

if diff_size $out_dir/08-affiliation_postprocessed.biom $expected_dir/08-affiliation_postprocessed.biom 0
then
	echo "difference in affiliation_postprocess :08-affiliation_postprocessed.biom" >&2
fi	

if diff_line $out_dir/08-affiliation_postprocessed.fasta $expected_dir/08-affiliation_postprocessed.fasta 0
then
	echo "difference in affiliation_postprocess :08-affiliation_postprocessed.fasta" >&2
fi	

echo "Step normalisation `date`"

if $run_programs
then
	normalisation.py \
	 -n 100 \
	 --input-biom $expected_dir/08-affiliation_postprocessed.biom \
	 --input-fasta $expected_dir/08-affiliation_postprocessed.fasta \
	 --output-biom $out_dir/09-normalisation.biom \
	 --output-fasta $out_dir/09-normalisation.fasta \
	 --summary $out_dir/09-normalisation.html \
	 --log-file $out_dir/09-normalisation.log
	 
	if [ $? -ne 0 ]
	then
		echo "Error in normalisation" >&2
		exit 1;
	fi
fi


##difficile à tester à cause du tirage aléatoire
if diff_line $out_dir/09-normalisation.fasta $expected_dir/09-normalisation.fasta 0
then
	echo "Difference in normalisation : 09-normalisation.fasta" >&2
fi

if diff_size $out_dir/09-normalisation.biom $expected_dir/09-normalisation.biom 0
then
	echo "Difference in normalisation : 09-normalisation.biom" >&2
fi

if diff_line $out_dir/09-normalisation.html $expected_dir/09-normalisation.html 0
then
	echo "Difference in normalisation : 09-normalisation.html" >&2
fi


echo "Step clusters_stat `date`"

if $run_programs
then
	clusters_stat.py \
	 --input-biom $expected_dir/09-normalisation.biom \
	 --output-file $out_dir/10-clustersStat.html \
	 --log-file $out_dir/10-clustersStat.log

	if [ $? -ne 0 ]
	then
		echo "Error in clusters_stat" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/10-clustersStat.html $expected_dir/10-clustersStat.html 0
then
	echo "Difference in clusters_stat : 10-clustersStat.html" >&2
fi

echo "Step affiliations_stat `date`"

if $run_programs
then
	affiliations_stat.py \
	 --input-biom $expected_dir/09-normalisation.biom \
	 --output-file $out_dir/11-affiliationsStat.html \
	 --log-file $out_dir/11-affiliationsStat.log \
	 --tax-consensus-tag "blast_taxonomy" \
	 --identity-tag "perc_identity" \
	 --coverage-tag "perc_query_coverage" \
	 --multiple-tag "blast_affiliations" \
	 --rarefaction-ranks Family Genus Species \
	 --taxonomic-ranks Domain Phylum Class Order Family Genus Species ## this is the default value of this option

	if [ $? -ne 0 ]
	then
		echo "Error in affiliations_stat" >&2
		exit 1;
	fi
fi

## tirage aléatoire sur les courbes de raréfaction
if diff_line $out_dir/11-affiliationsStat.html $expected_dir/11-affiliationsStat.html 1
then
	echo "Difference in nb diff affiliations_stat : 11-affiliationsStat.html" >&2
fi

if diff_size $out_dir/11-affiliationsStat.html $expected_dir/11-affiliationsStat.html 10
then
	echo "Difference in size affiliations_stat : 11-affiliationsStat.html" >&2
fi

echo "Step biom_to_tsv `date`"

if $run_programs
then
	biom_to_tsv.py \
	 --input-biom $expected_dir/09-normalisation.biom \
	 --input-fasta $expected_dir/09-normalisation.fasta \
	 --output-tsv $out_dir/12-biom2tsv.tsv \
	 --output-multi-affi $out_dir/12-biom2tsv-affiliation_multihit.tsv \
	 --log-file $out_dir/12-biom2tsv.log

	if [ $? -ne 0 ]
	then
		echo "Error in biom_to_tsv" >&2
		exit 1;
	fi
fi

if diff_line $out_dir/12-biom2tsv.tsv $expected_dir/12-biom2tsv.tsv 0
then
	echo "Difference in biom_to_tsv : 12-biom2tsv.tsv" >&2
fi

if diff_line $out_dir/12-biom2tsv-affiliation_multihit.tsv $expected_dir/12-biom2tsv-affiliation_multihit.tsv 0
then
	echo "Difference in biom_to_tsv : 12-affiliation_multihit.tsv" >&2
fi

echo "Step biom_to_stdBiom `date`"

if $run_programs
then
	biom_to_stdBiom.py \
	 --input-biom $expected_dir/09-normalisation.biom \
	 --output-biom $out_dir/13-affiliation_std.biom \
	 --output-metadata $out_dir/13-affiliation_multihit.tsv \
	 --log-file $out_dir/13-biom2stdbiom.log

	if [ $? -ne 0 ]
	then
		echo "Error in biom_to_stdBiom" >&2
		exit 1;
	fi
fi

if diff_size $out_dir/13-affiliation_std.biom $expected_dir/13-affiliation_std.biom 0
then
	echo "Difference in biom_to_stdBiom: 13-affiliation_std.biom" >&2
fi

if diff_line $out_dir/13-affiliation_multihit.tsv  $expected_dir/13-affiliation_multihit.tsv 0
then
	echo "Difference in biom_to_stdBiom: 13-affiliation_multihit.tsv " >&2
fi

echo "Step tsv_to_biom `date`"

if $run_programs
then
	tsv_to_biom.py \
	 --input-tsv $expected_dir/12-biom2tsv.tsv \
	 --input-multi-affi $expected_dir/12-biom2tsv-affiliation_multihit.tsv \
	 --output-biom $out_dir/14-tsv2biom.biom \
	 --output-fasta $out_dir/14-tsv2biom.fasta \
	 --log-file $out_dir/14-tsv2biom.log 

	if [ $? -ne 0 ]
	then
		echo "Error in tsv_to_biom" >&2
		exit 1;
	fi
fi

if diff_size $out_dir/14-tsv2biom.biom $expected_dir/14-tsv2biom.biom 0
then
	echo "Difference in tsv_to_biom : 14-tsv2biom.biom" >&2
fi

if diff_line $out_dir/14-tsv2biom.fasta $expected_dir/14-tsv2biom.fasta 0
then
	echo "Difference in tsv_to_biom : 14-tsv2biom.fasta" >&2
fi

echo "Step tree : mafft `date`"

if $run_programs
then
	tree.py \
	 --nb-cpus $nb_cpu \
	 --input-sequences $expected_dir/09-normalisation.fasta \
	 --biom-file $expected_dir/09-normalisation.biom \
	 --out-tree $out_dir/15-tree-mafft.nwk \
	 --html $out_dir/15-tree-mafft.html \
	 --log-file $out_dir/15-tree-mafft.log --debug

	if [ $? -ne 0 ]
	then
		echo "Error in tree : mafft" >&2
		exit 1;
	fi
fi

# mafft produit des résultats différents même avec la même version. Faire des tests individuel pour fasttree et phangorn
# mafft_aln=`ls $out_dir/*_mafft_aligned.fasta`
# pref=`echo $mafft_aln | sed 's/_mafft_aligned.fasta//'`
# FastTree -nt -gtr $mafft_aln > ${prefix}_fasttree.nwk2 2> ${prefix}_fasttree.stderr2
# ls -l ${prefix}_fasttree.nwk*
# sdiff -s ${prefix}_fasttree.nwk* |wc -l
# root_tree.R ${prefix}_fasttree.nwk2 $out_dir/15-tree-mafft.nwk2
# ls -l res_3.2.3-r4/15-tree-mafft.nwk*
# sdiff -s res_3.2.3-r4/15-tree-mafft.nwk* |wc -l


if diff_size $out_dir/15-tree-mafft.nwk $expected_dir/15-tree-mafft.nwk 10
then
	echo "Difference in tree : 15-tree-mafft.nwk" >&2
fi

if diff_line $out_dir/15-tree-mafft.html $expected_dir/15-tree-mafft.html 1
then
	echo "Difference in number of line in tree : 15-tree-mafft.html" >&2
fi

if diff_size $out_dir/15-tree-mafft.html $expected_dir/15-tree-mafft.html 10
then
	echo "Difference in in size in tree : 15-tree-mafft.html" >&2
fi

echo "Step phyloseq_import_data `date`"

if $run_programs
then
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
fi

if diff_size $out_dir/16-phylo_import.Rdata $expected_dir/16-phylo_import.Rdata 5
then
	echo "Difference in size phyloseq_import_data : 16-phylo_import.Rdata" >&2
fi


echo "Step phyloseq_composition `date`"

if $run_programs
then
	phyloseq_composition.py  \
	 --varExp EnvType --taxaRank1 Kingdom --taxaSet1 Bacteria --taxaRank2 Phylum --numberOfTaxa 9 \
	 --rdata $expected_dir/16-phylo_import.Rdata \
	 --html $out_dir/17-phylo_composition.nb.html \
	 --log-file $out_dir/17-phylo_composition.log

	 
	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_composition " >&2
		exit 1;
	fi
fi

grep Actinobacteria $out_dir/17-phylo_composition.nb.html | sed 's/},/},\n/g'  > tmp
grep Actinobacteria $expected_dir/17-phylo_composition.nb.html | sed 's/},/},\n/g'  > tmp2
if diff_line tmp tmp2 3
then
	echo "Difference in phyloseq_composition : 17-phylo_composition.nb.html" >&2
fi

rm tmp tmp2

echo "Step phyloseq_alpha_diversity `date`"

if $run_programs
then
	phyloseq_alpha_diversity.py  \
	 --varExp EnvType \
	 --rdata $expected_dir/16-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
	 --alpha-out $out_dir/18-phylo_alpha_div.tsv \
	 --html $out_dir/18-phylo_alpha_div.nb.html \
	 --log-file $out_dir/18-phylo_alpha_div.log

	 
	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_alpha_diversity " >&2
		exit 1;
	fi
fi

if diff_line $out_dir/18-phylo_alpha_div.tsv $expected_dir/18-phylo_alpha_div.tsv 0
then
	echo "Difference in phyloseq_alpha_diversity : 18-phylo_alpha_div.tsv" >&2
fi

echo "Step phyloseq_beta_diversity `date`"

if $run_programs
then
	phyloseq_beta_diversity.py  \
	 --varExp EnvType --distance-methods cc,unifrac \
	 --rdata $expected_dir/16-phylo_import.Rdata \
	 --matrix-outdir $out_dir \
	 --html $out_dir/19-phylo_beta_div.nb.html \
	 --log-file $out_dir/19-phylo_beta_div.log
	 
	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_beta_diversity " >&2
		exit 1;
	fi
fi

if diff_line $out_dir/unifrac.tsv $expected_dir/unifrac.tsv 0
then
	echo "Difference in phyloseq_beta_diversity : unifrac.tsv" >&2
fi

echo "Step phyloseq_structure `date`"

if $run_programs
then
	phyloseq_structure.py  \
	 --varExp EnvType --ordination-method MDS \
	 --rdata $expected_dir/16-phylo_import.Rdata --distance-matrix $expected_dir/unifrac.tsv \
	 --html $out_dir/20-phylo_structure.nb.html \
	 --log-file $out_dir/20-phylo_structure.log

	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_structure " >&2
		exit 1;
	fi
fi

a=`grep Axis.1 $out_dir/20-phylo_structure.nb.html |sed 's/Axis.1/Axis.1\n/'  | grep -v -m 1 'Axis.1' | sed 's/"/\n/' | head -n 1 |awk '{print $1}'`
b=`grep Axis.1 $expected_dir/20-phylo_structure.nb.html |sed 's/Axis.1/Axis.1\n/'  | grep -v -m 1 'Axis.1' | sed 's/"/\n/' | head -n 1 |awk '{print $1}'`
if [ $a != $b ]
then
	echo "Error in phylo_structure :  Diffrent Axis.1 contribution" >&2
fi

a=`grep Axis.2 $out_dir/20-phylo_structure.nb.html |sed 's/Axis.2/Axis.2\n/'  | grep -v -m 1 'Axis.2' | sed 's/"/\n/' | head -n 1 |awk '{print $1}'`
b=`grep Axis.2 $expected_dir/20-phylo_structure.nb.html |sed 's/Axis.2/Axis.2\n/'  | grep -v -m 1 'Axis.2' | sed 's/"/\n/' | head -n 1 |awk '{print $1}'`
if [ $a != $b ]
then
	echo "Difference in phylo_structure :  Diffrent Axis.2 contribution" >&2
fi

echo "Step phyloseq_clustering `date`"

if $run_programs
then
	phyloseq_clustering.py  \
	 --varExp EnvType \
	 --rdata $expected_dir/16-phylo_import.Rdata --distance-matrix $expected_dir/unifrac.tsv \
	 --html $out_dir/21-phylo_clustering.nb.html \
	 --log-file $out_dir/21-phylo_clustering.log

	 
	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_clustering " >&2
		exit 1;
	fi
fi

if diff_line $out_dir/21-phylo_clustering.nb.html $expected_dir/21-phylo_clustering.nb.html 0
then
	echo "Difference in phyloseq_clustering : 21-phylo_clustering.nb.html" >&2
fi

echo "Step phyloseq_manova `date`"

if $run_programs
then
	phyloseq_manova.py  \
	 --varExp EnvType \
	 --rdata $expected_dir/16-phylo_import.Rdata --distance-matrix $expected_dir/unifrac.tsv \
	 --html $out_dir/22-phylo_manova.nb.html \
	 --log-file $out_dir/22-phylo_manova.log

	if [ $? -ne 0 ]
	then
		echo "Error in phyloseq_manova " >&2
		exit 1;
	fi
fi

a=`grep EnvType $out_dir/22-phylo_manova.nb.html |grep -v adonis |sed 's/\*//g'`
b=`grep EnvType $expected_dir/22-phylo_manova.nb.html |grep -v adonis |sed 's/\*//g'`

if [ "$a" != "$b" ]
then
	echo "Difference in phyloseq_manova :  different statistical test result" >&2
fi

echo "Step deseq2_preprocess `date`"

if $run_programs
then
	deseq2_preprocess.py \
	 --data $expected_dir/16-phylo_import.Rdata \
	 --log-file $out_dir/23-deseq2_preprocess.log \
	 --out-Rdata $out_dir/23-deseq2_preprocess.Rdata \
	 --var EnvType

	if [ $? -ne 0 ]
	then
		echo "Error in deseq2_preprocess " >&2
		exit 1;
	fi
fi

if diff_line $out_dir/23-deseq2_preprocess.Rdata $expected_dir/23-deseq2_preprocess.Rdata 50
then
	echo "Difference in deseq2_preprocess : 23-deseq2_preprocess.Rdata " >&2
fi

echo "Step deseq2_visualisation `date`"

if $run_programs
then
	deseq2_visualisation.py \
	 --phyloseqData $expected_dir/16-phylo_import.Rdata \
	 --dds $expected_dir/23-deseq2_preprocess.Rdata \
	 --log-file $out_dir/24-deseq2_visualisation.log \
	 --html $out_dir/24-deseq2_visualisation.nb.html \
	 --var EnvType --mod1 BoeufHache --mod2 SaumonFume
	                            

	if [ $? -ne 0 ]
	then
		echo "Error in deseq2_visualisation " >&2
		exit 1;
	fi
fi

# récupérer les tableau CSV via la page HTML
# et faire sdiff

grep otu_01582 $out_dir/24-deseq2_visualisation.nb.html | sed 's/],/],\n/g' > /tmp/tmp
grep otu_01582 $expected_dir/24-deseq2_visualisation.nb.html | sed 's/],/],\n/g'  > /tmp/tmp1

if diff_line /tmp/tmp /tmp/tmp1 1
then
	echo "Difference in deseq2_visualisation : 24-deseq2_visualisation.nb.html  " >&2
fi

rm /tmp/tmp /tmp/tmp1
echo "Completed with success"

