#!/usr/bin/env python2.7
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
__version__ = 'r3.0-v1.0'
__email__ = 'frogs@inra.fr'
__status__ = 'dev'

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
def parse_blast(blast_file, query_list, ref_list):

    last_query = ""
    last_score = 0
    cpt = 0
    FH_in = open(blast_file)
    for line in FH_in:
        query = line.split("\t")[0].split("_FROGS")[0]
        #only for query to align
        if query not in query_list:
            continue
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

def extract_ref(input_fasta, input_blast_R1, input_blast_R2, input_ref, output_ref):

    # query_id list
    query_ids = []
    FH_in = FastaIO(input_fasta)
    for record in FH_in:
        query_ids.append(record.id.replace("_FROGS_combined",""))
    FH_in.close()

    # best ref id list
    best_ref = []

    # parse blast R1
    parse_blast(input_blast_R1, query_ids, best_ref)
    # parse blast R2
    parse_blast(input_blast_R2, query_ids, best_ref)

    # extract ref
    if len(best_ref) > 0 :
        FH_in = FastaIO(input_ref)
        FH_out = FastaIO(output_ref,"w")
        for record in FH_in:
            if record.id in best_ref:
                FH_out.write(record)
    
        FH_out.close()
        FH_in.close()

    return len(best_ref)

def submit_cmd( cmd, stdout_path, stderr_path):
    """
    @summary: Submits a command on system.
    @param cmd: [list] The command.
    @param stdout_path: The path to the file where the standard outputs will be written.
    @param stderr_path: The path to the file where the error outputs will be written.
    """
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # write down the stdout
    stdoh = open(stdout_path, "w")
    stdoh.write(stdout)
    stdoh.close()

    # write down the stderr
    stdeh = open(stderr_path, "w")
    stdeh.write(stderr)
    stdeh.close()

    # check error status
    if p.returncode != 0:
        stdeh = open(stderr_path)
        error_msg = "".join( map(str, stdeh.readlines()) )
        stdeh.close()
        raise StandardError( error_msg )

def get_needleall_version():
    """
    @summary: Return the emboss version.
    @return: [str] The emboss version.
    """
    version = None
    try:
        cmd = ["needlall", "-version"]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        version = stderr.strip()
    except:
        version = "unknown"
    return version

def process_needleall(ref_fasta, query_fasta, output_sam, needle_out, needle_err, log_file):

    cmd = ["needleall", "-asequence", ref_fasta, "-bsequence", query_fasta, "-outfile", output_sam, "-aformat3", "sam", "-gapopen", "10.0", "-gapextend", "0.5"]
    FH_log = Logger( log_file )
    FH_log.write("# Needlall version: " + get_needleall_version() + "\n")
    FH_log.write("# Needlall command: " + " ".join(cmd) + "\n")
    submit_cmd( cmd, needle_out, needle_err )
    FH_log.close()

def process(params):

    tmpFiles = TmpFiles( os.path.split(params.output_sam)[0] )
    extract_fasta_ref = tmpFiles.add(os.path.basename(params.reference))
    needle_out = tmpFiles.add(os.path.basename(params.output_sam)+"needle.out")
    needle_err = tmpFiles.add(os.path.basename(params.output_sam)+"needle.err")
    
    # Extract best blast ref
    try :
        Logger.static_write(params.log_file, "# Parsing blast alignment results to reduce reference databse\n")
        nb_ref = extract_ref(params.query_fasta, params.query_blast_R1, params.query_blast_R2, params.reference,extract_fasta_ref)

        # Launch needleall
        if nb_ref > 0 :
            Logger.static_write(params.log_file, "\tReducing reference databases to " + str(nb_ref) + " sequences\n\n")
            process_needleall(extract_fasta_ref, params.query_fasta, params.output_sam, needle_out, needle_err, params.log_file)
        else:
            Logger.static_write(params.log_file, "\tNeedle alignment not performed because no blast alignment available\n")
            open(params.output_sam, 'w').close()
    finally:
        if not params.debug:
            tmpFiles.deleteAll()

    
###################################################################################################################
###                                              MAIN                                                           ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Reduced reference fasta file based on best blast before runing NeedleAll alignement")
    parser.add_argument( '-v', '--version', action='version', version=__version__ + " [needleall " + get_needleall_version() + "]")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")

    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-r', '--reference', required=True, help='Reference fasta file')
    group_input.add_argument('-f', '--query-fasta', required=True, help='Query fasta file')
    group_input.add_argument('--query-blast-R1', required=True, help='Query R1 blast tsv file. NCBI Blastn+ need to be used with -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\' option)')
    group_input.add_argument('--query-blast-R2', required=True, help='Query R2 blast tsv file. NCBI Blastn+ need to be used with -outfmt \'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\' option)')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-s', '--output-sam', default="needleall.sam", help='Sam Needlall output file. [Default: %(default)s]')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)