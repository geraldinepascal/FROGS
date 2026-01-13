#!/bin/bash

nb_cpu=1
out_dir=res_check
expected_dir=res_dev
run_programs=true     ## if true lance les python sinon, fait uniquement les comparatifs de résultats

## Set ENV
## export PATH=../app:$PATH
## or
## conda activate frogsfunc@XXX
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

echo "Step frogsfunc_placeseqs `date`"

if $run_programs
then
	frogsfunc_placeseqs.py \
		--input-fasta data/frogsfunc.fasta \
		--input-biom data/frogsfunc.biom \
		--nb-cpus $nb_cpus \
		--placement-tool sepp \
		--output-tree  $out_dir/25-frogsfunc_placeseqs_tree.nwk \
		--excluded $out_dir/25-frogsfunc_placeseqs_excluded.tsv \
		--output-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
		--output-biom $out_dir/25-frogsfunc_placeseqs.biom \
		--closests-ref $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt \
		--log-file $out_dir/25-frogsfunc_placeseqs.log \
		--output-marker-copy $out_dir/25-frogsfunc_placeseqs_marker.tsv \
		--html $out_dir/25-frogsfunc_placeseqs_summary.html

	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_placeseqs " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_tree.nwk $expected_dir/25-frogsfunc_placeseqs_tree.nwk 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_tree.nwk" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_excluded.txt $expected_dir/25-frogsfunc_placeseqs_excluded.txt 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_excluded.txt" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs.fasta $expected_dir/25-frogsfunc_placeseqs.fasta 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs.fasta" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs.biom  $expected_dir/25-frogsfunc_placeseqs.biom 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs.biom" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt $expected_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_closests_ref_sequences.txt" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_marker.tsv $expected_dir/25-frogsfunc_placeseqs_marker.tsv 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_marker.txt" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_summary.html  $expected_dir/25-frogsfunc_placeseqs_summary.html 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_summary.html" >&2
fi

echo "Step frogsfunc_functions `date`"

if $run_programs
then
	frogsfunc_functions.py \
		--input-biom $out_dir/25-frogsfunc_placeseqs.biom \
		--input-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
		--input-tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
		--nb-cpus $nb_cpus \
		--marker-type 16S \
		--input-marker-copy $out_dir/25-frogsfunc_placeseqs_marker.tsv \
		--prefix-function-abund $out_dir/26-frogsfunc_functions_unstrat \
		--output-asv-copy-norm $out_dir/26-frogsfunc_functions_marker_norm.tsv \
		--output-weighted $out_dir/26-frogsfunc_functions_weighted_nsti.tsv \
		--output-excluded $out_dir/26-frogsfunc_functions_excluded.tsv \
		--output-fasta $out_dir/26-frogsfunc_functions.fasta \
		--output-biom $out_dir/26-frogsfunc_functions.biom \
		--log-file $out_dir/26-frogsfunc_functions.log \
		--html $out_dir/26-frogsfunc_functions_summary.html

	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_functions " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/26-frogsfunc_functions_unstrat_EC.tsv $expected_dir/26-frogsfunc_functions_unstrat_EC.tsv 0
then
	echo "Difference in frogsfunc_functions : 26-frogsfunc_functions_unstrat_EC.tsv " >&2
fi

if diff_line $out_dir/EC_copynumbers_predicted.tsv $expected_dir/EC_copynumbers_predicted.tsv 0
then
	echo "Difference in frogsfunc_functions : EC_copynumbers_predicted.tsv " >&2
fi


if diff_line $out_dir/26-frogsfunc_functions_marker_norm.tsv $expected_dir/26-frogsfunc_functions_marker_norm.tsv 0
then
	echo "Difference in frogsfunc_functions : 26-frogsfunc_functions_marker_norm.tsv " >&2
fi

if diff_line $out_dir/26-frogsfunc_functions_weighted_nsti.tsv $expected_dir/26-frogsfunc_functions_weighted_nsti.tsv 0
then
	echo "Difference in frogsfunc_functions : 26-frogsfunc_functions_weighted_nsti.tsv " >&2
fi

if diff_line $out_dir/26-frogsfunc_functions_excluded.tsv $expected_dir/26-frogsfunc_functions_excluded.tsv 0
then
	echo "Difference in frogsfunc_functions : 26-frogsfunc_functions_excluded.tsv " >&2
fi

if diff_line $out_dir/26-frogsfunc_functions_summary.html $expected_dir/26-frogsfunc_functions_summary.html 0
then
	echo "Difference in frogsfunc_functions : 26-frogsfunc_functions_summary.html " >&2
fi


echo "Step frogsfunc_pathways `date`"

if $run_programs
then	
	frogsfunc_pathways.py \
		--input-tsv $out_dir/26-frogsfunc_functions_unstrat_EC.tsv \
		--nb-cpus $nb_cpus \
		--normalisation \
		--strat-contrib \
		--input-asv-copy-norm $out_dir/26-frogsfunc_functions_marker_norm.tsv \
		--input-fun-copy $out_dir/EC_copynumbers_predicted.tsv \
		--output-pathways-abund $out_dir/27-strat-frogsfunc_pathways_unstrat.tsv \
		--output-pathways-contrib $out_dir/27-strat-frogsfunc_pathways_strat.tsv \
		--output-pathways-predictions $out_dir/27-strat-frogsfunc_pathways_predictions.tsv \
		--output-pathways-abund-per-seq $out_dir/27-strat-frogsfunc_pathways_unstrat_per_seq.tsv \
		--log-file $out_dir/27-strat-frogsfunc_pathways.log \
		--html  $out_dir/27-strat-frogsfunc_pathways_summary.html

	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_pathways " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/27-strat-frogsfunc_pathways_unstrat.tsv $expected_dir/27-strat-frogsfunc_pathways_unstrat.tsv 0
then
	echo "Difference in frogsfunc_pathways : 27-strat-frogsfunc_pathways_unstrat.tsv " >&2
fi

if diff_line $out_dir/27-strat-frogsfunc_pathways_strat.tsv $expected_dir/27-strat-frogsfunc_pathways_strat.tsv 0
then
	echo "Difference in frogsfunc_pathways : 27-strat-frogsfunc_pathways_strat.tsv " >&2
fi


if diff_line $out_dir/27-strat-frogsfunc_pathways_predictions.tsv $expected_dir/27-strat-frogsfunc_pathways_predictions.tsv 0
then
	echo "Difference in frogsfunc_pathways : 27-strat-frogsfunc_pathways_predictions.tsv " >&2
fi


if diff_line $out_dir/27-strat-frogsfunc_pathways_unstrat_per_seq.tsv $expected_dir/27-strat-frogsfunc_pathways_unstrat_per_seq.tsv 0
then
	echo "Difference in frogsfunc_pathways : 27-strat-frogsfunc_pathways_unstrat_per_seq.tsv " >&2
fi


if diff_line $out_dir/27-strat-frogsfunc_pathways_summary.html $expected_dir/27-strat-frogsfunc_pathways_summary.html 0
then
	echo "Difference in frogsfunc_pathways : 27-strat-frogsfunc_pathways_summary.html " >&2
fi

echo "Completed with success"
