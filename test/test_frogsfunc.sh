#!/bin/bash
out_dir=$1

# Check parameters
if [ "$#" -ne 1 ]; then
    echo "ERROR: Illegal number of parameters." ;
    echo 'Command usage: test_frogsfunc.sh <OUT_FOLDER>' ;
    exit 1 ;
fi

# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi

echo "Step frogsfunc_placeseqs `date`"

frogsfunc_placeseqs.py \
 --input-fasta data/frogsfunc.fasta \
 --input-biom data/frogsfunc.biom \
 --placement-tool sepp \
 --output-tree  $out_dir/25-frogsfunc_placeseqs_tree.nwk \
 --excluded $out_dir/25-frogsfunc_placeseqs_excluded.txt \
 --output-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
 --output-biom $out_dir/25-frogsfunc_placeseqs.biom \
 --closests-ref $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt \
 --log-file $out_dir/25-frogsfunc_placeseqs.log \
 --output-marker $out_dir/25-frogsfunc_marker.tsv \
 --summary $out_dir/25-frogsfunc_placeseqs_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_placeseqs " >&2
    exit 1;
fi

echo "Step frogsfunc_functions `date`"

frogsfunc_functions.py \
 --input-biom $out_dir/25-frogsfunc_placeseqs.biom \
 --input-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
 --input-tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
 --marker-type 16S \
 --input-marker $out_dir/25-frogsfunc_marker.tsv \
 --output-function-abund $out_dir/26-frogsfunc_functions_unstrat.tsv \
 --output-otu-norm $out_dir/26-frogsfunc_functions_marker_norm.tsv \
 --output-weighted $out_dir/26-frogsfunc_functions_weighted_nsti.tsv \
 --output-excluded $out_dir/26-frogsfunc_functions_excluded.txt \
 --output-fasta $out_dir/26-frogsfunc_function.fasta \
 --output-biom $out_dir/26-frogsfunc_function.biom \
 --log-file $out_dir/26-frogsfunc_functions.log \
 --summary $out_dir/26-frogsfunc_functions_summary.html


if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_functions " >&2
    exit 1;
fi

echo "Step frogsfunc_pathways `date`"

frogsfunc_pathways.py \
 --input-file $out_dir/26-frogsfunc_functions_unstrat_EC.tsv \
 --normalisation \
 --per-sequence-contrib \
 --per-sequence-abun $out_dir/26-frogsfunc_functions_marker_norm.tsv \
 --per-sequence-function $out_dir/EC_copynumbers_predicted.tsv \
 --output-pathways-abund $out_dir/27-frogsfunc_pathways_unstrat.tsv \
 --output-pathways-contrib $out_dir/27-frogsfunc_pathways_strat.tsv \
 --output-pathways-predictions $out_dir/27-frogsfunc_pathways_predictions.tsv \
 --output-pathways-abund-per-seq $out_dir/27-frogsfunc_pathways_unstrat_per_seq.tsv \
 --log-file $out_dir/27-frogsfunc_pathways.log \
 --summary  $out_dir/27-frogsfunc_pathways_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_pathways " >&2
    exit 1;
fi

echo "Completed with success"