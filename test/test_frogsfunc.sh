#!/bin/bash
frogs_dir=$1
out_dir=$2

# Check parameters
if [ "$#" -ne 2 ]; then
    echo "ERROR: Illegal number of parameters." ;
    echo 'Command usage: test_frogsfunc.sh <FROGS_FOLDER> <OUT_FOLDER>' ;
    exit 1 ;
fi

# Set ENV
export PATH=$frogs_dir/app:$PATH

# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi

echo "Step frogsfunc_placeseqs `date`"

frogsfunc_placeseqs.py \
 --input-fasta $frogs_dir/test/data/frogsfunc.fasta \
 --input-biom $frogs_dir/test/data/frogsfunc.biom \
 --placement-tool sepp \
 --out-tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
 --excluded $out_dir/25-frogsfunc_placeseqs_excluded.txt \
 --insert-fasta $out_dir/25-frogsfunc_placeseqs.fasta \
 --insert-biom $out_dir/25-frogsfunc_placeseqs.biom \
 --closests-ref $out_dir/25-frogsfunc_placeseqs_closests_ref_sequences.txt \
 --log-file $out_dir/25-frogsfunc_placeseqs.log \
 --html $out_dir/25-frogsfunc_placeseqs_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_placeseqs " >&2
    exit 1;
fi

echo "Step frogsfunc_copynumbers `date`"

frogsfunc_copynumbers.py \
 --input-biom $out_dir/25-frogsfunc_placeseqs.biom \
 --tree $out_dir/25-frogsfunc_placeseqs_tree.nwk \
 --output-marker $out_dir/26-frogsfunc_copynumbers_marker.tsv \
 --output-function $out_dir/26-frogsfunc_copynumbers_predicted_functions.tsv \
 --log-file $out_dir/26-frogsfunc_copynumbers.log \
 --html $out_dir/26-frogsfunc_copynumbers_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_copynumbers " >&2
    exit 1;
fi

echo "Step frogsfunc_functions `date`"

frogsfunc_functions.py \
 --input-biom $out_dir/25-frogsfunc_placeseqs.biom \
 --function $out_dir/26-frogsfunc_copynumbers_predicted_functions.tsv \
 --marker $out_dir/26-frogsfunc_copynumbers_marker.tsv \
 --function-abund $out_dir/27-frogsfunc_functions_unstrat.tsv \
 --seqtab $out_dir/27-frogsfunc_functions_marker_norm.tsv \
 --weighted $out_dir/27-frogsfunc_functions_weighted_nsti.tsv \
 --excluded $out_dir/27-frogsfunc_functions_excluded.txt \
 --log-file $out_dir/27-frogsfunc_functions.log \
 --html $out_dir/27-frogsfunc_functions_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_functions " >&2
    exit 1;
fi

echo "Step frogsfunc_pathways `date`"

frogsfunc_pathways.py \
 --input-file $out_dir/27-frogsfunc_functions_unstrat.tsv \
 --pathways-abund $out_dir/28-frogsfunc_pathways_unstrat.tsv \
 --log-file $out_dir/28-frogsfunc_pathways.log \
 --html $out_dir/28-frogsfunc_pathways_summary.html --debug

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_pathways " >&2
    exit 1;
fi

echo "Completed with success"
