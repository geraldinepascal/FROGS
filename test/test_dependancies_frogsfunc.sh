#!/bin/bash

out_dir=res_4.0.1_to_check
expected_dir=frogsfunc_4.0.0
run_programs=false     ## if true lance les python sinon, fait uniquement les comparatifs de résultats

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
	 --placement-tool sepp \
	 --out-tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
	 --excluded $out_dir/25-frogsfunc_placeseqs_excluded.txt \
	 --insert-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
	 --insert-biom $out_dir/25-frogsfunc_placeseqs.biom \
	 --closests-ref $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt \
	 --html $out_dir/25-frogsfunc_placeseqs_report.html \
	 --log-file $out_dir/25-frogsfunc_placeseqs.log

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

if diff_size $out_dir/25-frogsfunc_placeseqs.biom  $expected_dir/25-frogsfunc_placeseqs.biom 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs.biom" >&2
fi

if diff_line $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt $expected_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_closests_ref_sequences.txt" >&2
fi

if diff_size $out_dir/25-frogsfunc_placeseqs_report.html  $expected_dir/25-frogsfunc_placeseqs_report.html 0
then
	echo "difference in frogsfunc_placeseqs, mask mode :25-frogsfunc_placeseqs_report.html" >&2
fi

echo "Step frogsfunc_copynumbers `date`"

if $run_programs
then
	frogsfunc_copynumbers.py \
	 --input-biom $out_dir/25-frogsfunc_placeseqs.biom \
	 --tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
	 --output-marker $out_dir/26-frogsfunc_copynumbers_marker.tsv \
	 --output-function $out_dir/26-frogsfunc_copynumbers_predicted_functions.tsv \
	 --html $out_dir/26-frogsfunc_copynumbers_report.html \
	 --log-file $out_dir/26-frogsfunc_copynumbers.log
	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_copynumbers " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/26-frogsfunc_copynumbers_marker.tsv $expected_dir/26-frogsfunc_copynumbers_marker.tsv 0
then
	echo "Difference in frogsfunc_copynumbers : 26-frogsfunc_copynumbers_marker.tsv " >&2
fi

if diff_line $out_dir/26-frogsfunc_copynumbers_predicted_functions.tsv $expected_dir/26-frogsfunc_copynumbers_predicted_functions.tsv 0
then
	echo "Difference in frogsfunc_copynumbers : 26-frogsfunc_copynumbers_predicted_functions.tsv " >&2
fi

if diff_line $out_dir/26-frogsfunc_copynumbers_report.html $expected_dir/26-frogsfunc_copynumbers_report.html 0
then
	echo "Difference in frogsfunc_copynumbers : 26-frogsfunc_copynumbers_report.html " >&2
fi


echo "Step frogsfunc_functions `date`"

if $run_programs
then
	frogsfunc_functions.py \
	 --input-biom $out_dir/25-frogsfunc_placeseqs.biom \
	 --function $out_dir/26-frogsfunc_copynumbers_predicted_functions.tsv \
	 --marker $out_dir/26-frogsfunc_copynumbers_marker.tsv \
	 --function-abund $out_dir/27-frogsfunc_functions_unstrat.tsv \
	 --seqtab $out_dir/27-frogsfunc_functions_marker_norm.tsv \
	 --weighted $out_dir/27-frogsfunc_functions_weighted_nsti.tsv \
	 --excluded $out_dir/27-frogsfunc_functions_excluded.txt \
	 --html $out_dir/27-frogsfunc_functions_report.html \
	 --log-file $out_dir/27-frogsfunc_functions.log

	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_genefamilies " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/27-frogsfunc_functions_unstrat.tsv $expected_dir/27-frogsfunc_functions_unstrat.tsv 0
then
	echo "Difference in frogsfunc_functions : 27-frogsfunc_functions_unstrat.tsv " >&2
fi

if diff_line $out_dir/27-frogsfunc_functions_marker_norm.tsv $expected_dir/27-frogsfunc_functions_marker_norm.tsv 0
then
	echo "Difference in frogsfunc_functions : 27-frogsfunc_functions_marker_norm.tsv " >&2
fi

if diff_line $out_dir/27-frogsfunc_functions_weighted_nsti.tsv $expected_dir/27-frogsfunc_functions_weighted_nsti.tsv 0
then
	echo "Difference in frogsfunc_functions : 27-frogsfunc_functions_weighted_nsti.tsv " >&2
fi

if diff_line $out_dir/27-frogsfunc_functions_excluded.txt $expected_dir/27-frogsfunc_functions_excluded.txt 0
then
	echo "Difference in frogsfunc_functions : 27-frogsfunc_functions_excluded.txt " >&2
fi

if diff_line $out_dir/27-frogsfunc_functions_report.html $expected_dir/27-frogsfunc_functions_report.html 0
then
	echo "Difference in frogsfunc_functions : 27-frogsfunc_functions_report.html " >&2
fi


echo "Step frogsfunc_pathways `date`"

if $run_programs
then
	frogsfunc_pathways.py \
	 --input-file $out_dir/27-frogsfunc_functions_unstrat.tsv \
	 --pathways-abund $out_dir/28-frogsfunc_pathways_unstrat.tsv \
	 --html $out_dir/28-frogsfunc_pathways_report.html \
	 --log-file $out_dir/28-frogsfunc_pathways.log

	if [ $? -ne 0 ]
	then
	    echo "Error in frogsfunc_pathways " >&2
	    exit 1;
	fi
fi

if diff_line $out_dir/28-frogsfunc_pathways_unstrat.tsv $expected_dir/28-frogsfunc_pathways_unstrat.tsv 0
then
	echo "Difference in frogsfunc_pathways : 28-frogsfunc_pathways_unstrat.tsv " >&2
fi

if diff_line $out_dir/28-frogsfunc_pathways_report.html $expected_dir/28-frogsfunc_pathways_report.html 0
then
	echo "Difference in frogsfunc_pathways : 28-frogsfunc_pathways_report.html " >&2
fi

echo "Completed with success"
