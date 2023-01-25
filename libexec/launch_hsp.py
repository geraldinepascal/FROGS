#!/usr/bin/env python3
#
# Copyright (C) 2022 INRAE
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

__author__ = ' Vincent Darbot - GENPHYSE '
__copyright__ = 'Copyright (C) 2022 INRAE'
__license__ = 'GNU General Public License'
__version__ = '4.0.1'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'dev'

import os
import sys
import pandas as pd
import argparse
from subprocess import Popen, PIPE
import threading
import multiprocessing
from multiprocessing import Queue

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPAT
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
# Default table PATH

from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class HspMarker(Cmd):
    """
    @summary: Predict number of marker copies (16S, 18S or ITS) for each cluster sequence (i.e OTU).
    """
    def __init__(self, observed_marker_table, in_tree, hsp_method, output, log):
        """
        @param observed_marker_table: [str] Path to marker table file if marker studied is not 16S.
        @param in_tree: [str] Path to resulting tree file with inserted clusters sequences from frogsfunc_placeseqs.
        @param hsp_method: [str] HSP method to use.
        @param output: [str] PICRUSt2 marker output file.
        """
        if observed_marker_table is None:
            input_marker = " -i 16S"
        else:
            input_marker = " --observed_trait_table " + observed_marker_table        

        Cmd.__init__(self,
                 'hsp.py',
                 'predict marker copy number per ASV sequence.', 
                 input_marker + " -t " + in_tree + " --hsp_method " + hsp_method + " -o " + output + " --calculate_NSTI  2> " + log,
                "--version")

        self.output = output

    def get_version(self):
        return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

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
    stdoh = open(stdout_path, "wt")
    stdoh.write(stdout.decode('utf-8'))
    stdoh.close()

    # write down the stderr
    stdeh = open(stderr_path, "wt")
    stdeh.write(stderr.decode('utf-8'))
    stdeh.close()

    # check error status
    if p.returncode != 0:
        stdeh = open(stderr_path,'rt')
        error_msg = "".join( map(str, stdeh.readlines()) )
        stdeh.close()
        raise_exception( Exception( "\n\n#ERROR : " + error_msg + "\n\n" ))


def process_hsp_function(in_traits, observed_trait_table, in_tree, hsp_method, outputs, logs):
    if type(in_traits) != list:
        in_traits = [in_traits]
    if type(outputs) != list:
        outputs = [outputs]
    if type(logs) != list:
        logs = [logs]
    # run hsp.py
    for idx, in_trait in enumerate(in_traits):
        if observed_trait_table is None:
            input_function = " --in_trait " + in_trait
            message = "## Process function : " + in_trait + "\n"
        else:
            input_function = " --observed_trait_table " + observed_trait_table
            message = "## Process function table : " + observed_trait_table + "\n"
        FH_log = Logger( logs[idx] )
        FH_log.write(message)
        cmd = ["hsp.py", input_function.split()[0], input_function.split()[1] ,"-t", in_tree, "--hsp_method", hsp_method, "-o", outputs[idx]]
        FH_log.write("## hsp.py command: " + " ".join(cmd) + "\n")
        submit_cmd( cmd, logs[idx], logs[idx] )
        # FH_log.write("".join(open(tmp_out).readlines()) + "\n" )
        # FH_log.write("".join(open(tmp_err).readlines()) + "\n" )
        FH_log.close()


def parallel_submission( function, inputs, tree, hsp_method, outputs, logs, cpu_used):
    processes = [{'process':None, 'inputs':None, 'tree':tree, 'hsp_method':hsp_method, 'outputs':None, 'log_files':None} for trait in range(cpu_used)]
    # Launch processes
    for trait in range(len(inputs)):
        process_idx = trait % cpu_used
        processes[process_idx]['inputs'] = inputs[trait]
        processes[process_idx]['outputs'] = outputs[trait]
        processes[process_idx]['log_files'] = logs[trait]
    
    for current_process in processes:
        if trait == 0:  # First process is threaded with parent job
            
            current_process['process'] = threading.Thread(target=function,
                                                          args=(current_process['inputs'], None, tree, hsp_method, current_process['outputs'], current_process['log_files']))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                          args=(current_process['inputs'], None, tree, hsp_method, current_process['outputs'], current_process['log_files']))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise_exception( Exception("\n\n#ERROR : Error in sub-process execution.\n\n"))

def append_results(functions_outputs, marker_file, logs_hsp, final_output, log_file):
    """
    """
    otu_to_nsti = dict()
    otu_to_results = dict()
    with open(marker_file) as FH_marker:
        FH_marker.readline()
        for li in FH_marker:
            li = li.strip().split('\t')
            otu_to_nsti[li[0]] = li[2]
            otu_to_results[li[0]] = list()

    for current_file in functions_outputs:
        with open(current_file) as FH_input:
            i_functions = dict()
            header = FH_input.readline().strip().split('\t')[1:]
            for i in range(len(header)):
                i_functions[i] = header[i]
            
            for li in FH_input:
                otu = li.strip().split('\t')[0]
                abundances = li.strip().split('\t')[1:]
                for cur_abund in range(len(abundances)):
                    otu_to_results[otu].append({ i_functions[cur_abund] : str(abundances[cur_abund])})
    
    FH_out = open(final_output , "wt")
    for otu, functions in otu_to_results.items():
        FH_out.write("sequence\t" + "\t".join(list(function.keys())[0] for function in functions) + "\n")
        break
    for otu, functions in otu_to_results.items():
        FH_out.write(otu + "\t" + "\t".join(list(function.values())[0] for function in functions) + "\n")
    FH_out.close()
    # Append log
    FH_log = Logger(log_file)
    FH_log.write("\n")
    for current_file in logs_hsp:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write(line)
        FH_input.close()
        FH_log.write("\n")
    FH_log.write("\n")
    FH_log.close()


def task_marker(args):
    # Check for ITS or 18S input
    if args.marker_type in ["ITS", "18S"]:
        if args.input_marker_table is None:
            parser.error("\n\n#ERROR : --input-marker-table required when studied marker is not 16S!\n\n")
    args.to_launch = "marker"
    return args


def task_function(args):
    # Check for 16S input
    if args.marker_type == "16S" and (not 'EC' in args.functions and not 'KO' in args.functions):
            parser.error("\n\n#ERROR : --functions : 'EC' and/or 'KO' must be at least indicated (others functions are optionnal)")
    # Check for ITS or 18S input
    if args.marker_type in ["ITS", "18S"]:
        if args.input_function_table is None:
            parser.error("\n\n#ERROR : --input-function-table required when studied marker is not 16S!\n\n")
        else:
            args.functions = None
    args.to_launch = "function"
    return args

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='predict marker of gene copy number' )
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    subparsers = parser.add_subparsers()
    # Inputs
    parser_marker = subparsers.add_parser('marker', help='Predict marker copy number per ASV sequence.')
    parser_marker.add_argument('-m', '--marker-type', required=True, choices=['16S','ITS','18S'], help='Marker gene to be analyzed.')
    parser_marker.add_argument('-t', '--input-tree', required=True, type=str, help='frogsfunc_placeseqs output tree in newick format containing both studied sequences (i.e. ASVs or OTUs) and reference sequences.')
    parser_marker.add_argument('--hsp-method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use. mp: predict discrete traits using max parsimony. emp_prob: predict discrete traits based on empirical state probabilities across tips. subtree_average: predict continuous traits using subtree averaging. pic: predict continuous traits with phylogentic independent contrast. scp: reconstruct continuous traits using squared-change parsimony (default: %(default)s).')
    group_input_marker_other = parser_marker.add_argument_group( 'ITS and 18S ' )
    group_input_marker_other.add_argument('--input-marker-table',help="The input marker table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. (ex $PICRUSt2_PATH/default_files/fungi/ITS_counts.txt.gz). Required.")
    group_output_marker = parser_marker.add_argument_group( 'Outputs' )
    group_output_marker.add_argument('-o', '--output-marker', default="frogsfunc_copynumbers_marker.tsv", type=str, help='Output table of predicted marker gene copy numbers per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.[Default: %(default)s]')
    parser_marker.set_defaults(func=task_marker)

    parser_function = subparsers.add_parser('function', help='Predict gene copy number per ASV sequence.')
    parser_function.add_argument('-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    parser_function.add_argument('-m', '--marker-type', required=True, choices=['16S','ITS','18S'], help='Marker gene to be analyzed.')
    parser_function.add_argument('-i', '--marker-file', required=True, type=str, help='Table of predicted marker gene copy numbers (frogsfunc_placeseqs output : frogsfunc_marker.tsv).')
    parser_function.add_argument('-t', '--input-tree', required=True, type=str, help='frogsfunc_placeseqs output tree in newick format containing both studied sequences (i.e. ASVs or OTUs) and reference sequences.')
    parser_function.add_argument('--hsp-method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use. mp: predict discrete traits using max parsimony. emp_prob: predict discrete traits based on empirical state probabilities across tips. subtree_average: predict continuous traits using subtree averaging. pic: predict continuous traits with phylogentic independent contrast. scp: reconstruct continuous traits using squared-change parsimony (default: %(default)s).')
    group_input_16S = parser_function.add_argument_group( '16S' )
    group_input_16S.add_argument('-f', '--functions', default=["EC"], nargs='+', choices=['EC', 'KO', 'COG', 'PFAM', 'TIGRFAM','PHENO'], help="Specifies which function databases should be used (%(default)s). EC is used by default because necessary for frogsfunc_pathways. At least EC or KO is required. To run the command with several functions, separate the functions with spaces (ex: -i EC PFAM).")
    group_input_other = parser_function.add_argument_group( 'ITS and 18S ' )
    group_input_other.add_argument('--input-function-table',help="The path to input functions table describing directly observed functions, in tab-delimited format.(ex $PICRUSt2_PATH/default_files/fungi/ec_ITS_counts.txt.gz). Required.")
    group_output_function = parser_function.add_argument_group( 'Outputs' )
    group_output_function.add_argument('--output-dir', default="frogsfunc_function_results", help='Output directory for function predictions.')
    group_output_function.add_argument('-o', '--output-function', default="frogsfunc_copynumbers_functions.tsv", type=str, help='Output table with predicted function abundances per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.[Default: %(default)s]')
    parser_function.set_defaults(func=task_function)

    # Output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='List of commands executed.')
    args = parser.parse_args()
    args.to_launch = str()
    args.func(args)
    prevent_shell_injections(args)
    tmp_files=TmpFiles(os.path.split(args.input_tree)[0])

    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        if args.to_launch == "marker":
            tmp_hsp_marker = tmp_files.add( 'tmp_hsp_marker.log' )
            HspMarker(args.input_marker_table, args.input_tree, args.hsp_method, args.output_marker, tmp_hsp_marker).submit(args.log_file)

        if args.to_launch == "function":
            tmp_hsp_function = tmp_files.add( 'tmp_hsp_function.log' )

            # if args.functions is not None:
            suffix_name = "_copynumbers_predicted.tsv"
            functions_outputs = [args.output_dir + "/" + trait + "_copynumbers_predicted.tsv" for trait in args.functions]
            logs_hsp = [tmp_files.add( trait + "_tmp_hsp_function.log") for trait in args.functions]
            if len(args.functions) == 1 or args.nb_cpus == 1:
                Logger.static_write(args.log_file, '\n\nRunning ' + " ".join(args.functions) + ' functions prediction.\n')
                process_hsp_function(args.functions, args.input_function_table, args.input_tree, args.hsp_method, functions_outputs, logs_hsp)

            else:
                parallel_submission( process_hsp_function, args.functions, args.input_tree, args.hsp_method, functions_outputs, logs_hsp, len(args.functions) )
            
            # else:
            #     tmp_output_function = tmp_files.add( "copynumbers_predicted.tsv")
            #     process_hsp_function(args.functions, args.input_function_table, args.input_tree, args.hsp_method, args.output_function, tmp_hsp_function)
    
    finally:
        if not args.debug:
            tmp_files.deleteAll()
