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
__version__ = '1.0.0'
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
        print "\t[FROGS CMD]:\t" + cmd
        subprocess.check_call( cmd, shell=True )
    else:
        print "\t[FROGS CMD]:\t" + cmd + " > " + output
        subprocess.check_call( cmd + " > " + output, shell=True )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    print "[SOFTWARE]:\tFROGS"
    start_time = time.time()

    # Manage parameters
    parser = argparse.ArgumentParser(description='Launch FROGS workflow.')
    parser.add_argument( '--min-length', type=int, required=True, help='The minimum amplicon length.')
    parser.add_argument( '--max-length', type=int, required=True, help='The maximum amplicon length.')
    parser.add_argument( '--already-contiged', action='store_true', default=False, help='The archive contains 1 file by sample : Reads 1 and Reads 2 are already contiged by pair.' )
    parser.add_argument( '--without-primers', action='store_true', default=False, help='Use this option when you use custom sequencing primers and these primers are the PCR primers. In this case the reads do not contain the PCR primers.' )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help='The maximum number of CPUs used. [Default: 1]')
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-a', '--affiliation-databank', required=True, help='The path to the affiliation databank (format: fasta).')
    group_input.add_argument('-r', '--reads-archive', help='The path to the input archive (format: tar).')
    group_input.add_argument('-s', '--samples-files', nargs='*', help='The path to samples sequences files (format: fastq).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-b', '--output-biom', required=True, help='The final abundance table with affiliation (format: BIOM).')
    group_output.add_argument('-f', '--output-fasta', required=True, help='The final sequences (format: fasta).')
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output_biom)

    preprocess_options = ""
    if args.without_primers: preprocess_options += " --without-primers"
    if args.already_contiged: preprocess_options += " --already-contiged"
    input_options = ""
    if args.reads_archive is not None:
        input_options = " --input-archive " + args.reads_archive
    elif args.samples_files is not None and len(args.samples_files) != 0:
        input_options = " --input-R1 " + " ".join(args.samples_files)
    else:
        raise Exception("Option '--reads-archive' or '--samples-files' must be provided.")
    exec_cmd( "preprocess.py illumina" + \
    " --nb-cpus " + str(args.nb_cpus) + \
    " --min-amplicon-size " + str(args.min_length) + \
    " --max-amplicon-size " + str(args.max_length) + \
    preprocess_options + \
    input_options + \
    " --output-dereplicated " + os.path.join(out_dir, "prepro.fasta") + \
    " --output-count " + os.path.join(out_dir, "prepro.tsv") + \
    " --summary " + os.path.join(out_dir, "prepro_summary.html") + \
    " --log-file " + os.path.join(out_dir, "prepro_log.txt") )


    exec_cmd( "clustering.py" + \
    " --nb-cpus " + str(args.nb_cpus) + \
    " --input-fasta " + os.path.join(out_dir, "prepro.fasta") + \
    " --input-count " + os.path.join(out_dir, "prepro.tsv") + \
    " --output-biom " + os.path.join(out_dir, "clustering.biom") + \
    " --output-fasta " + os.path.join(out_dir, "clustering.fasta") + \
    " --output-compo " + os.path.join(out_dir, "clustering_compo.tsv") + \
    " --log-file " + os.path.join(out_dir, "clustering_log.txt") + \
    " --distance 3" + \
    " --denoising" )


    exec_cmd( "remove_chimera.py" + \
    " --nb-cpus " + str(args.nb_cpus) + \
    " --input-fasta " + os.path.join(out_dir, "clustering.fasta") + \
    " --input-biom " + os.path.join(out_dir, "clustering.biom") + \
    " --non-chimera " + os.path.join(out_dir, "removeChimera.fasta") + \
    " --out-abundance " + os.path.join(out_dir, "removeChimera.biom") + \
    " --summary " + os.path.join(out_dir, "removeChimera_summary.html") + \
    " --log-file " + os.path.join(out_dir, "removeChimera_log.txt") )


    exec_cmd( "filters.py" + \
    " --input-biom " + os.path.join(out_dir, "removeChimera.biom") + \
    " --input-fasta " + os.path.join(out_dir, "removeChimera.fasta") + \
    " --output-fasta " + args.output_fasta + \
    " --output-biom " + os.path.join(out_dir, "filters.biom") + \
    " --excluded " + os.path.join(out_dir, "filters_excluded.txt") + \
    " --summary " + os.path.join(out_dir, "filters_summary.html") + \
    " --log-file " + os.path.join(out_dir, "filters_log.txt") + \
    " --min-abundance 0.00005" )


    exec_cmd( "affiliation_OTU.py" + \
    " --nb-cpus " + str(args.nb_cpus) + \
    " --reference " + args.affiliation_databank + \
    " --input-fasta " + args.output_fasta + \
    " --input-biom " + os.path.join(out_dir, "filters.biom") + \
    " --output-biom " + args.output_biom + \
    " --summary " + os.path.join(out_dir, "affiliationOTU_summary.html") + \
    " --log-file " + os.path.join(out_dir, "affiliationOTU_log.txt") + \
    " --java-mem 20" )


    end_time = time.time()
    print "\t[FROGS EXEC_TIME]:\t" + str(end_time - start_time)
