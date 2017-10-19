#!/usr/bin/env python2.7
#
# Copyright (C) 2014 INRA
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

__author__ = 'Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs@inra.fr'
__status__ = 'prod'

import os
import sys
import time
import argparse
import subprocess
from subprocess import Popen, PIPE
import threading
import multiprocessing
from multiprocessing import Queue

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsSequenceIO import *
from frogsUtils import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def get_ITSx_version():
    """
    @summary: Return the ITSx version.
    @return: [str] The ITSx version.
    """
    version = None
    try:
        cmd = ["ITSx", "-h"]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        for line in stderr.split("\n"):
            if line.startswith("Version"):
                return line.split()[-1].strip()
    except:
        version = "unknown"
    return version

def get_fasta_nb_seq( fasta_file ):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(fasta_file):
        FH_input = open( fasta_file )
    else:
        FH_input = gzip.open( fasta_file )
    nb_seq = 0
    for line in FH_input:
        if line.startswith(">"):
            nb_seq += 1
    FH_input.close()
    return nb_seq

def split_fasta(fasta_file, tmp_files_manager, nb_file, out_list, log_file):
    """
    @summary: split fasta in nb_file and returne outfile list
    @param fasta_file: [str] Path to the fasta file to process.
    @param tmp_files_manager: [TmpFiles] The temporary file manager.
    @param nb_file : [int] The number of file to generate.
    @param out_list : [list] List of output fasta file
    @param log_file : [srt] path to logfile
    """
    out_files = list()
    record_iter = FastaIO(fasta_file)
    for idx, record in enumerate(record_iter):
        out_file_idx = idx % nb_file
        if len(out_files) == 0 or not out_file_idx < len(out_files):
            new_out_file = tmp_files_manager.add( os.path.basename(fasta_file) + "_" + str(out_file_idx) )
            out_files.append({ 'file_path': new_out_file,
                               'file_handle': FastaIO(new_out_file, "w"),
                               'nb_seq': 0
            })
            out_list.append( os.path.abspath(new_out_file) )
        out_files[out_file_idx]['nb_seq'] += 1
        out_files[out_file_idx]['file_handle'].write( record )
    for out_file in out_files:
        out_file['file_handle'].close()

    # Log
    FH_log = Logger(log_file)
    FH_log.write("# split " + fasta_file + " in " + str(nb_file) + " fasta files\n")
    FH_log.write("Results\n")
    for out_file in out_files:
        FH_log.write( "\tWrote " + str(out_file['nb_seq']) + " records to " + out_file['file_path'] + "\n" )
    FH_log.close()

def submit_cmd( cmd, cwd=None):
    """
    @summary: Submits a command on system.
    @param cmd: [list] The command.
    @param stdout_path: The path to the file where the standard outputs will be written.
    @param stderr_path: The path to the file where the error outputs will be written.
    """
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=cwd)
    stdout, stderr = p.communicate()

    # check error status
    if p.returncode != 0:
        # stdeh = open(stderr)
        error_msg = "".join( map(str, stderr.readlines()) )
        # stdeh.close()
        raise StandardError( error_msg )

def parallel_submission( function, inputs, cwds, outputs, logs, cpu_used):
    processes = [{'process':None, 'inputs':None, 'cwd' : None, 'outputs':None, 'log_files':None} for idx in range(cpu_used)]
    # Launch processes
    for idx in range(len(inputs)):
        process_idx = idx % cpu_used
        processes[process_idx]['inputs'] = inputs[idx]
        processes[process_idx]['cwd'] = cwds[idx]
        processes[process_idx]['outputs'] = outputs[idx]
        processes[process_idx]['log_files'] = logs[idx]

    for current_process in processes:
        if idx == 0:  # First process is threaded with parent job
            current_process['process'] = threading.Thread(target=function,
                                                          args=(current_process['inputs'], current_process['cwd'], current_process['outputs'], current_process['log_files']))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                                 args=(current_process['inputs'], current_process['cwd'], current_process['outputs'], current_process['log_files']))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise Exception("Error in sub-process execution.")

def replaceTag(input, tag1, tag2, out, log):
    """
    @summary : modify concatenation tag from ITSx and remove 5.8S first or last base
    @param input  : [str] Path to ITS1 sequence fasta file
    @param tag1     : [str] the sequence tag to replace
    @param tag2     : [str] the sequence tag to use
    @param out     : [str] Path to fasta concatenated ITS sequence fasta file
    """
    FH_in = FastaIO(input)
    FH_out = FastaIO(out,"w")
    ITS1=  0
    ITS2 = 0
    ITS1_2= 0 
    for record in FH_in:

        splitted = record.string.split(tag1)
        if len(splitted[0]) == 1 :
            record.string = splitted[1]
            if not ";size=" in record.id : 
                record.id += "_ITS2"
            else :
                record.id = record.id.replace(";size=", "_ITS2;size=")
            ITS2 += 1
        elif len(splitted[1]) == 1:
            record.string = splitted[0]
            if not ";size=" in record.id : 
                record.id += "_ITS1"
            else :
                record.id = record.id.replace(";size=", "_ITS1;size=")
            ITS1 += 1
        else:
            record.string = record.string.replace(tag1,tag2)
            if not "FROGS_combined" in record.id :
                if not ";size=" in record.id : 
                    record.id += "_FROGS_combined_ITS1_ITS2"
                else :
                    record.id = record.id.replace(";size=", "_FROGS_combined_ITS1_ITS2;size=")
            else:
                if not ";size=" in record.id : 
                    record.id += "_ITS1_ITS2"
                else :
                    record.id = record.id.replace(";size=", "_ITS1_ITS2;size=")
            ITS1_2 += 1 

        FH_out.write(record)
    FH_in.close()
    FH_out.close()
    FH_log = Logger( log )
    FH_log.write("## replace ITSx \"-----\" tag combining ITS 1 and 2 by 100 N and remove firts or last 5.8S base when only one ITS region is included\n")
    FH_log.write("##Results\n")
    FH_log.write("Output file :" + os.path.split(out)[1] + "\n" + 
                     "nb_ITS1_only: " + str(ITS1)  + "\n" +
                     "nb_ITS2_only: " + str(ITS2)  + "\n" +
                     "nb_ITS1_ITS2: " + str(ITS1_2)  + "\n" )
    FH_log.close()

def process_ITSx(in_fasta, cwd, out, log_file):

    os.mkdir(cwd)
    prefix = os.path.splitext(os.path.split(in_fasta)[1])[0]

    # run ITSx
    FH_log = Logger( log_file )
    FH_log.write("## Input file : " + os.path.split(in_fasta)[1] + "\n" ) 
    FH_log.write("## in working directory: " + cwd + "\n")
    cmd = ["ITSx", "-i", in_fasta, "-o", prefix , "--preserve", "T","--concat","T","-t","F"]
    FH_log.write("## ITSx command: " + " ".join(cmd) + "\n")
    submit_cmd( cmd , cwd )
    FH_log.close()

    replaceTag(os.path.join(cwd,prefix+".concat.fasta"), "-----", 100*"N", out, log_file)

def append_results(appended_fasta, appended_log, fasta_out, log_file):
    """
    @summary: Summeriza log files in one and concatenated multiprocessed fasta file.
    @param appended_fasta: [list] List of fasta files to append.
    @param appended_log: [list] List of log files to summarized.
    @param fasta_out: [str] The fasta file where contents of others fasta are appended.
    @param log_file: [str] The log file where contents of others log are summed.
    """
    # Append fasta
    FH_fasta = FastaIO(fasta_out , "w")
    for current_file in appended_fasta:
        FH_input = FastaIO(current_file)
        for record in FH_input:
            FH_fasta.write(record)
        FH_input.close()
    FH_fasta.close()

    # Append log
    FH_log = Logger(log_file)
    FH_log.write("\n")
    for current_file in appended_log:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write(line)
        FH_input.close()
        FH_log.write("\n")
    FH_log.write("\n")
    FH_log.close()


def Update_count(input_count, input_fasta, output_count) :
    """
    @summary : update count by keeping only line from the resulting fasta file
    @param input_count [str] : Path to input count file
    @param input_fasta [str] : Path to ITSx filtered fasta file
    @param output_count [str] : Path to updated count file
    """
    
    seq_dict = dict()

    FH_in = FastaIO(input_fasta)
    for record in FH_in:
        seq_dict[record.id.split(";size=")[0].replace("_ITS2","").replace("_ITS1","")] = record.id
    FH_in.close()

    FH_in = open(input_count)
    FH_out = open(output_count, "w")

    for line in FH_in : 
        if line.startswith("#") :
            FH_out.write(line)
        else:
            seq_id = line.split()[0]
            if seq_id in seq_dict :
                FH_out.write( seq_dict[seq_id]+ "\t" + "\t".join(line.split()[1:]) + "\n" )
            elif seq_id + "_FROGS_combined" in seq_dict :
                FH_out.write( seq_dict[seq_id + "_FROGS_combined"] + "\t" + "\t".join(line.split()[1:]) + "\n")

def main_process(args):

    try : 
        
        tmpFiles = TmpFiles(os.path.split(args.output_fasta)[0])

        # set number of process
        nb_seq = get_fasta_nb_seq(args.input_fasta)
        
        if args.nb_cpus == 1 or nb_seq < 10:
            in_fasta = os.path.abspath(args.input_fasta)
            tmp_dir = tmpFiles.add_dir(os.path.split(args.output_fasta)[1])
            process_ITSx(in_fasta, tmp_dir, args.output_fasta, args.log_file)

        else:
            fasta_ITSx_list = list()
            ITSx_outputs = list()
            logs_ITSx = list()
            tmp_dirs = list()

            split_fasta(args.input_fasta, tmpFiles, max(1,args.nb_cpus-1), fasta_ITSx_list, args.log_file)
            ITSx_outputs = [tmpFiles.add(os.path.basename(current_fasta) + ".fasta") for current_fasta in fasta_ITSx_list]
            logs_ITSx = [tmpFiles.add(os.path.basename(current_fasta) + "_itsx.log") for current_fasta in fasta_ITSx_list]
            tmp_dirs = [ tmpFiles.add_dir(os.path.split(current_fasta)[1]) for current_fasta in fasta_ITSx_list ]
            parallel_submission( process_ITSx, fasta_ITSx_list, tmp_dirs, ITSx_outputs, logs_ITSx, len(fasta_ITSx_list) )

            # Logs
            append_results(ITSx_outputs, logs_ITSx, args.output_fasta, args.log_file)

        Update_count(args.input_count, args.output_fasta, args.output_count)
    
    finally:
        if not args.debug:
            tmpFiles.deleteAll()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Search for ITS region per subset of sequences'
    )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ + " [ITSx " + get_ITSx_version() + "]" )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-fasta', required=True, help='The fasta input sequences to treat' )
    group_input.add_argument( '-c', '--input-count', required=True, help='The count tsv file associated with input fasta file' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-fasta', default='ITS.fasta', help='Fasta with ITS sequences only. [Default: %(default)s]' )
    group_output.add_argument( '-a', '--output-count', default='ITS.count.tsv', help='Updated count tsv file. [Default: %(default)s]' )
    group_output.add_argument( '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.' )
    args = parser.parse_args()

    # Process
    main_process(args)