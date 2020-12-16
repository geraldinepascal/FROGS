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

__author__ = 'Sigenae INRA Jouy en Josas'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import time
import argparse
import subprocess
from frogsSequenceIO import *

def exec_cmd( cmd, output=None ):
    """
    @summary: Execute command line and check status.
    @param cmd: [str] Command line.
    @param output: [str] Path to the file where the stdout will be written.
    """
    if output is None:
        print cmd
        subprocess.check_call( cmd, shell=True )
    else:
        print cmd + " > " + output
        subprocess.check_call( cmd + " > " + output, shell=True )

if __name__ == "__main__":

    start_time = time.time()
    print "[SOFTWARE]:\tQiime\tstart time : ",start_time
    # Manage parameters
    parser = argparse.ArgumentParser(description="Launch qiime workflow.")
    parser.add_argument( '-n', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: 1]")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-folder', required=True, help='The path to the input fastq files folder.')
    group_input.add_argument('-r', '--ref-fasta', required=True, help='The path to the reference fasta file')
    group_input.add_argument('-t', '--ref-tax', help='The path to the reference tax file')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', required=True, help='The final abundance table with affiliation (format: BIOM).')
    group_output.add_argument('--output-fasta', required=True, help='The final sequences (format: fasta).')
    args = parser.parse_args()

    working_path_prefix = os.path.join(os.path.split(args.output_biom)[0], "qiime_workdir")
    
    # Launch qiime
    qiime_fasta=""
    qiime_biom=""
    
    # concatenate and filter reads
    file_list=",".join([os.path.join(args.input_folder,f) for f in os.listdir(args.input_folder) ] )
    sample_name_list=",".join([f.split("-")[0] for f in os.listdir(args.input_folder)])
    
    exec_cmd("qiime; split_libraries_fastq.py -i " + file_list \
        + " --sample_ids " + sample_name_list \
        +" -o " + os.path.join(working_path_prefix, "qiime_preprocess") \
        +" --barcode_type 'not-barcoded' " \
        +" --phred_offset 33")

    merge_fasta=os.path.join(working_path_prefix,"qiime_preprocess","seqs.fna")
    
    # Launch chimera identification (in Qiime)
    exec_cmd("qiime; identify_chimeric_seqs.py -i "+ merge_fasta \
        + " -m usearch61 --suppress_usearch61_ref " \
        + " -o " + os.path.join(working_path_prefix,"usearch61_chimeras") )
        
    # Remove chimera
    qiime_input_fasta=os.path.join(working_path_prefix,"usearch61_chimeras","seqs_chimeras_filtered.fna")
    exec_cmd("qiime; filter_fasta.py -f " + merge_fasta \
        + " -o " + qiime_input_fasta \
        + " -s " + os.path.join(working_path_prefix,"usearch61_chimeras","chimeras.txt") \
        + " -n")
        
    # Launch Qiime
    cpus_opt = ""
    if args.nb_cpus > 1 :
        cpus_opt = " -aO "+str(args.nb_cpus)
             
    qiime_command ="qiime; pick_open_reference_otus.py -i " + qiime_input_fasta\
        + cpus_opt \
        + " -o "+ os.path.join(working_path_prefix, "pick_open_reference_otus") \
        + " -r "+ args.ref_fasta \
        + " --suppress_align_and_tree" \
        + " --suppress_taxonomy_assignment"
    exec_cmd(qiime_command)
    
    qiime_fasta=os.path.join(working_path_prefix, "pick_open_reference_otus","rep_set.fna")
    if args.ref_tax is not None:
        exec_cmd("qiime; assign_taxonomy.py -o " + os.path.join(working_path_prefix,"uclust_assigned_taxonomy") \
            + " -i " + qiime_fasta \
            + " -t " + args.ref_tax \
            + " -r " + args.ref_fasta )

        exec_cmd("qiime_completTax.py -i " + os.path.join(working_path_prefix,"uclust_assigned_taxonomy","rep_set_tax_assignments.txt") \
            + " -o " + os.path.join(working_path_prefix,"uclust_assigned_taxonomy","rep_set_completeTax_assignments.txt") )

        qiime_biom=os.path.join(working_path_prefix, "otu_table_mc2_w_tax.biom")
        exec_cmd("biom add-metadata -i " + os.path.join(working_path_prefix, "pick_open_reference_otus","otu_table_mc2.biom") \
         + " --observation-metadata-fp " + os.path.join(working_path_prefix,"uclust_assigned_taxonomy","rep_set_completeTax_assignments.txt") \
         + " -o " + qiime_biom \
         + " --sc-separated taxonomy --observation-header OTUID,taxonomy ")
        
    else:
         qiime_biom=os.path.join(working_path_prefix, "pick_open_reference_otus","otu_table_mc2.biom")

        
    # convert to properly biom file
    exec_cmd("biom convert -i "+ qiime_biom \
        +" -o "+ args.output_biom \
        +" --table-type=\"OTU table\" --to-json" )
        
    # add true seed ref
    exec_cmd("addSeedsRef.py -s "+ qiime_fasta \
        + " -r " + file_list.replace(","," ") \
        + " -a " + args.output_fasta  )

    end_time = time.time()
    print "\t[QIIME EXEC_TIME]:\t" + str(end_time - start_time)

