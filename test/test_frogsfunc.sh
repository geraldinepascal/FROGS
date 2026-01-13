#!/bin/bash
out_dir=$1
nb_cpus=$2

# conda activate frogsfunc@5.0.2
# export PATH=/home/maria/workspace/FROGS_dev/app:$PATH


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
    --nb-cpus $nb_cpus \
    --placement-tool sepp \
    --output-tree  $out_dir/25-frogsfunc_placeseqs_tree.nwk \
    --excluded $out_dir/25-frogsfunc_placeseqs_excluded.txt \
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

echo "Step frogsfunc_functions `date`"

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

echo "Step frogsfunc_pathways `date`"
    # --strat-contrib need additionnal input --input-asv-copy-norm --input-fun-copy
frogsfunc_pathways.py \
    --input-tsv $out_dir/26-frogsfunc_functions_unstrat_EC.tsv \
    --nb-cpus $nb_cpus \
    --output-pathways-abund $out_dir/27-frogsfunc_pathways_unstrat.tsv \
    --log-file $out_dir/27-frogsfunc_pathways.log \
    --html  $out_dir/27-frogsfunc_pathways_summary.html

frogsfunc_pathways.py \
    --input-tsv $out_dir/26-frogsfunc_functions_unstrat_EC.tsv \
    --nb-cpus $nb_cpus \
    --normalisation \
    --strat-contrib \
    --input-asv-copy-norm $out_dir/26-frogsfunc_functions_marker_norm.tsv \
    --input-fun-copy $out_dir/EC_copynumbers_predicted.tsv \
    --output-pathways-abund $out_dir/27-strat-norm-frogsfunc_pathways_unstrat.tsv \
    --output-pathways-contrib $out_dir/27-strat-norm-frogsfunc_pathways_strat.tsv \
    --output-pathways-predictions $out_dir/27-strat-norm-frogsfunc_pathways_predictions.tsv \
    --output-pathways-abund-per-seq $out_dir/27-strat-norm-frogsfunc_pathways_unstrat_per_seq.tsv \
    --log-file $out_dir/27-strat-norm-frogsfunc_pathways.log \
    --html  $out_dir/27-strat-norm-frogsfunc_pathways_summary.html

if [ $? -ne 0 ]
then
    echo "Error in frogsfunc_pathways " >&2
    exit 1;
fi

echo "Completed with success"



frogsfunc_pathways.py --nb-cpus ${GALAXY_SLOTS:-1} \
    --input-tsv '/tmp/tmpxdlpq188/files/f/9/1/dataset_f915b7f0-164b-4097-aee4-c896bd610cbf.dat' \
    --map /home/maria/workspace/FROGS-wrappers_dev/tools/frogs_func/test-data/databases/frogs_picrust2_db/pathway_mapfiles/metacyc_path2rxn_struc_filt_pro.txt   \
    --html /tmp/tmpxdlpq188/job_working_directory/000/2/outputs/dataset_0bb243c0-4ab1-4605-aa22-83e38875f101.dat \
    --output-pathways-abund '/tmp/tmpxdlpq188/job_working_directory/000/2/outputs/dataset_57486ae8-bc4a-497b-a350-1ee140083574.dat']
