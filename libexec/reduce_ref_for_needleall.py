#!/usr/bin/env python3
#
# Copyright (C) 2018 INRA
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

__author__ = 'Maria Bernard INRA - SIGENAE'
__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os, sys
import argparse
import re
import subprocess
from subprocess import Popen, PIPE

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
os.environ['PATH'] = CURRENT_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsSequenceIO import *

###################################################################################################################
###                                                 FUNCTION                                                    ###
###################################################################################################################
def parse_blast(blast_file, ref_list):

    last_query = ""
    last_score = 0
    cpt = 0
    FH_in = open(blast_file)
    for line in FH_in:
        query = line.split("\t")[0].split("_FROGS")[0]
        # init last query
        if last_query == "":
            last_query = query

        # for each query keep all references of the 100 best score
        while query == last_query  and line is not None and cpt <=100:
            ref = line.split("\t")[1]
            score = line.split("\t")[11]
            #count nb best scores
            if score != last_score:
                cpt += 1
                last_score = score
            # add ref 
            if ref not in ref_list:
                ref_list.append(ref)
            #next line
            line = next(FH_in, None)
            if line is not None:
                query = line.split("\t")[0].split("_FROGS")[0]

        #new query
        if line is not None:
            last_query = query
            last_score = line.split("\t")[11]
            ref = line.split("\t")[1]
            cpt = 1
            if ref not in ref_list:
                ref_list.append(ref)

    FH_in.close()

def extract_ref(input_blast_R1, input_blast_R2, input_ref, output_ref):

    # best ref id list
    best_ref = []

    # parse blast R1
    parse_blast(input_blast_R1, best_ref)
    # parse blast R2
    parse_blast(input_blast_R2, best_ref)

    # extract ref
    c = 0
    if len(best_ref) > 0 :
        FH_in = FastaIO(input_ref)
        FH_out = FastaIO(output_ref,"wt")
        for record in FH_in:
            c += 1
            if record.id in best_ref:
                FH_out.write(record)
    
        FH_out.close()
        FH_in.close()

    return c,len(best_ref)


def process(params):

    # Extract best blast ref

    Logger.static_write(params.log_file, "# Parsing blast alignment results to reduce reference databse\n")
    nb_tot, nb_ref = extract_ref(params.query_blast_R1, params.query_blast_R2, params.reference, params.output_fasta)
    Logger.static_write(params.log_file, "\tReducing reference databases from " + str(nb_tot) + " to " + str(nb_ref) + " sequences\n\n")

    
###################################################################################################################
###                                              MAIN                                                           ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Reduced reference fasta file based on R1 and R2 100 best score blast alignment")
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")

    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-r', '--reference', required=True, help='Reference fasta file')
    group_input.add_argument('--query-blast-R1', required=True, help='Query R1 blast tsv file. NCBI Blastn+ need to be used with -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\' option)')
    group_input.add_argument('--query-blast-R2', required=True, help='Query R2 blast tsv file. NCBI Blastn+ need to be used with -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\' option)')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-f', '--output-fasta', default="reduced_ref.fasta", help='Reduced reference fasta file [Default: %(default)s]')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)
