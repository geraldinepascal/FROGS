#!/usr/bin/env python2.7
#
# Copyright (C) 2016 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.1'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'beta'


import os
import sys
import time
import argparse
import subprocess


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def exec_cmd( cmd, output=None ):
    """
    @summary: Execute command line and check status.
    @param cmd: [str] Command line.
    @param output: [str] Path to the file where the stdout will be written.
    """
    if output is None:
        print "\t[Uparse CMD]:\t" + cmd
        subprocess.check_call( cmd, shell=True )
    else:
        print "\t[Uparse CMD]:\t" + cmd + " > " + output
        subprocess.check_call( cmd + " > " + output, shell=True )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    print "[SOFTWARE]:\tUparse"
    start_time = time.time()

    # Manage parameters
    parser = argparse.ArgumentParser(description="Launch uparse workflow.")
    parser.add_argument( '--min-length', type=int, required=True, help="The minimum amplicon length.")
    parser.add_argument( '--max-length', type=int, required=True, help="The maximum amplicon length.")
    parser.add_argument( '--already-contiged', action='store_true', default=False, help='Reads 1 and Reads 2 are already contiged by pair.')
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: 1]")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-d', '--databank', help='The path to the affiliation databank (format: udb).')
    group_input.add_argument('-f', '--input-folder', required=True, help='The path to the input folder.')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', required=True, help='The final abundance table with affiliation (format: BIOM).')
    group_output.add_argument('--output-fasta', required=True, help='The final sequences (format: fasta).')
    args = parser.parse_args()

    working_path_prefix = os.path.join(os.path.split(args.output_biom)[0], "uparse_tmp")
    reads = list()


    # Merge pairs
    if not args.already_contiged:
        exec_cmd( "usearch -fastq_mergepairs " \
            + os.path.join(args.input_folder, "*.fastq.gz") + " " \
            + "-fastqout " + working_path_prefix + "_merged.fastq " \
            + "-fastq_minmergelen " + str(args.min_length) + " " \
            + "-fastq_maxmergelen " + str(args.max_length) + " " \
            + "-threads " + str(args.nb_cpus) + " " \
            + "-relabel @" )
    else:
        labeled_fastq = list()
        idx = 0
        for sample_filename in os.listdir(args.input_folder):
            if sample_filename.endswith(".fastq"):
                idx += 1
                reads.append(os.path.join(args.input_folder, sample_filename))
                labeled_fastq.append(working_path_prefix  + "_spl_" + str(idx) + ".fastq")
                exec_cmd( "usearch -fastq_filter " \
                    + os.path.join(args.input_folder, sample_filename) + " " \
                    + "-fastqout " + working_path_prefix  + "_spl_" + str(idx) + ".fastq " \
                    + "-threads " + str(args.nb_cpus) + " " \
                    + "-relabel @" )
        exec_cmd( "cat " + " ".join(labeled_fastq) + " > " + working_path_prefix + "_merged.fastq " )


    # Filter sequences       
    exec_cmd( "usearch -fastq_filter " \
        + working_path_prefix + "_merged.fastq " \
        + "-fastaout " + working_path_prefix + "_filtered.fasta " \
        + "-fastq_maxee 1.0 " \
        + "-fastq_maxns 0 " \
        + "-threads " + str(args.nb_cpus) )


    # Dereplicate
    exec_cmd( "usearch -derep_fulllength " \
        + working_path_prefix + "_filtered.fasta " \
        + "-sizeout " \
        + "-fastaout " + working_path_prefix + "_uniques.fasta " \
        + "-threads " + str(args.nb_cpus) + " " )


    # Sort sequences and remove singletons
    exec_cmd( "usearch -sortbysize " \
        + working_path_prefix + "_uniques.fasta " \
        + "-fastaout " + working_path_prefix + "_sorted " \
        + "-minsize 2" )


    # Clustering
    exec_cmd( "usearch -cluster_otus " \
        + working_path_prefix + "_sorted " \
        + "-otus " + working_path_prefix + "_seeds.fasta " \
        + "-uparseout " + working_path_prefix + "_clusters.txt " \
        + "-relabel Cluster_ " \
        + "-sizein " \
        + "-sizeout" )


    # Remove chimera with ref
    #print "usearch -uchime_ref $1_seeds.fasta -db gold.fa -strand plus -nonchimeras $1_seeds_refChim.fasta -threads 1"
    #cmd = "usearch -uchime_ref $1_seeds.fasta -db gold.fa -strand plus -nonchimeras $1_seeds_refChim.fasta -threads 1"


    # Affiliation
    fasta_before_abund = working_path_prefix + "_seeds.fasta"
    if args.databank is not None:
        fasta_before_abund = working_path_prefix + "_affiliation.fasta"
        exec_cmd( "usearch -utax " \
            + working_path_prefix + "_seeds.fasta " \
            + "-db " + args.databank + " " \
            + "-fastaout " + fasta_before_abund + " " \
            + "-strand both " \
            + "-threads " + str(args.nb_cpus) )


    # Create abundance table
    """
     Problem 0.97
     Problem in BIOM creation when fasta contains affiliations: one OTU loose his ID
         exec_cmd( "usearch -usearch_global " \
             + working_path_prefix + "_merged.fastq " \
             + "-db " + fasta_before_abund + " " \
             + "-biomout " + args.output_biom + " " \
             + "-strand both " \
             + "-id 0.97 "  \
             + "-threads " + str(args.nb_cpus) )
    """
    exec_cmd( "usearch -usearch_global " \
        + working_path_prefix + "_merged.fastq " \
        + "-db " + working_path_prefix + "_seeds.fasta " \
        + "-biomout " + working_path_prefix + "_woAffi.biom " \
        + "-strand both " \
        + "-id 0.97 "  \
        + "-threads " + str(args.nb_cpus) )
    if args.databank is not None:
        exec_cmd( "addUtaxFromFasta.py " \
            + "--input-fasta " + fasta_before_abund + " "  \
            + "--input-biom " + working_path_prefix + "_woAffi.biom " \
            + "--output-biom " + args.output_biom + " " \
            + "--taxonomy-tag taxonomy" )
    else:
        exec_cmd( "ln -sf " + working_path_prefix + "_woAffi.biom " + args.output_biom )


    # Add reference ID in seeds descriptions and remove size
    ##gp##exec_cmd( "addSeedsRef.py " \
        ##gp##+ "--seeds-fasta " + working_path_prefix + "_seeds.fasta "  \
        ##gp##+ "--reads " + " ".join(reads) + " " \
        ##gp##+ "--annotated-fasta " + args.output_fasta )


    end_time = time.time()
    print "\t[Uparse EXEC_TIME]:\t" + str(end_time - start_time)
