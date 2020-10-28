#!/usr/bin/env python3
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
__email__ = 'frogs-support@inrae.fr'
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
from frogsBiom import BiomIO

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
                               'file_handle': FastaIO(new_out_file, "wt"),
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
        error_msg = "".join( map(str, stderr.decode('utf-8').readlines()) )
        # stdeh.close()
        raise Exception( "\n\n#ERROR : " + error_msg + "\n\n" )

def parallel_submission( function, inputs, its, cwds, outputs, logs, cpu_used):
    processes = [{'process':None, 'inputs':None, 'its':its, 'cwd' : None, 'outputs':None, 'log_files':None} for idx in range(cpu_used)]
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
                                                          args=(current_process['inputs'],current_process['its'], current_process['cwd'], current_process['outputs'], current_process['log_files']))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                                 args=(current_process['inputs'], current_process['its'], current_process['cwd'], current_process['outputs'], current_process['log_files']))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise Exception("\n\n#ERROR : Error in sub-process execution.\n\n")

def parseITSxResult(input_dir, prefix, its, out, log):
    """
    @summary : rename sequence of targeted ITS and summarize ITSx detections
    @param input  : [str] Path to ITSx output directory
    @param prefix     : [str] ITSx output files prefix 
    @param its     : [str] ITS region targeted (and to keep), either ITS1 or ITS2
    @param out     : [str] Path to fasta targeted ITS region sequence fasta file
    @param log     : [str] Path to log file to summarize ITSx detections
    """
    fasta_files = [ os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith(".fasta") ]
    count_ITSx = dict()

    FH_out = FastaIO(out,"wt")
    for fasta in fasta_files:
        # based on file name recover detection type : ITS1 ITS2 full chimeric no_detection.
        detection_type = os.path.basename(fasta).replace(prefix,"").replace(".fasta","")[1:] # first char is either "." or "_"
        count_ITSx[detection_type] = 0
        FH_in = FastaIO(fasta)
        for record in FH_in:
            count_ITSx[detection_type] += 1
            if detection_type == its :
                """
                if not ";size=" in record.id : 
                    record.id += "_"+its
                else :
                    record.id = record.id.replace(";size=", "_"+its+";size=")
                """
                FH_out.write(record)

        FH_in.close()

    FH_out.close()

    FH_log = Logger( log )
    FH_log.write("##Results\n")
    #FH_log.write("Output file :" + os.path.split(out)[1] + "\n" )
    for detection_type in count_ITSx :
        if its == 'no_detections':
            if detection_type == its:
                FH_log.write("\tnb "+detection_type+ " (removed): " + str(count_ITSx[detection_type]) + "\n")
            else : 
                FH_log.write("\tnb "+detection_type+ " (kept): " + str(count_ITSx[detection_type]) + "\n")
            
        else:
            if detection_type == its:
                FH_log.write("\tnb "+detection_type+ " (kept): " + str(count_ITSx[detection_type]) + "\n")
            else : 
                FH_log.write("\tnb "+detection_type+ " (removed): " + str(count_ITSx[detection_type]) + "\n")
    FH_log.close()

def process_ITSx(in_fasta, its, cwd, out, log_file):

    os.mkdir(cwd)
    prefix = os.path.splitext(os.path.split(in_fasta)[1])[0]

    # run ITSx
    FH_log = Logger( log_file )
    FH_log.write("## Input file : " + os.path.split(in_fasta)[1] + "\n" ) 
    FH_log.write("## in working directory: " + cwd + "\n")
    cmd = ["ITSx", "-i", in_fasta, "-o", prefix , "--preserve", "T","-t","F","--save_regions","all"]
    FH_log.write("## ITSx command: " + " ".join(cmd) + "\n")
    submit_cmd( cmd , cwd )
    FH_log.close()

    parseITSxResult(cwd, prefix, its, out, log_file)

def append_results(appended_fasta, appended_log, fasta_out, log_file):
    """
    @summary: Summeriza log files in one and concatenated multiprocessed fasta file.
    @param appended_fasta: [list] List of fasta files to append.
    @param appended_log: [list] List of log files to summarized.
    @param fasta_out: [str] The fasta file where contents of others fasta are appended.
    @param log_file: [str] The log file where contents of others log are summed.
    """
    # Append fasta
    FH_fasta = FastaIO(fasta_out , "wt")
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


def remove_observations( removed_observations, input_biom, output_biom ):
    """
    @summary: Removes the specified list of observations.
    @param removed_observations: [list] The names of the observations to remove.
    @param input_biom: [str] The path to the input BIOM.
    @param output_biom: [str] The path to the output BIOM.
    """
    biom = BiomIO.from_json( input_biom )
    biom.remove_observations( removed_observations )
    BiomIO.write( output_biom, biom )

def write_summary( samples_names, log_remove_global, log_remove_spl, out_file ):
    """
    @summary: Writes the summary file.
    @param samples_names: [list] The samples names.
    @param sample_logs: [list] list of sample logs files
    @param log_remove_global: [dict] The global remove metrics.
    @param log_remove_spl: [dict] The remove metrics by sample.
    @param out_file: [str] Path to the summary file.
    """
    
    """
    # Collect metrics
    detection_results = dict()
    for idx, sample in enumerate(samples_names):
        detection_results[sample] = get_sample_resuts(sample_logs[idx])
    """
    # Writes output
    FH_out = open(out_file, "wt")
    global_remove_results = [ log_remove_global['nb_removed'], log_remove_global['nb_kept'],
                              log_remove_global['abundance_removed'], log_remove_global['abundance_kept']]
    FH_out.write( '##Metrics global\n' )
    FH_out.write( "\t".join(['#Nb removed', 'Nb kept', 'Abundance removed', 'Abundance kept']) + "\n" )
    FH_out.write( "\t".join(map(str, global_remove_results)) + "\n" )
    FH_out.write( "\n" )

    FH_out.write( '##Metrics by sample\n' )
    FH_out.write( "\t".join(['#Sample name', 'Kept nb', 'Kept abundance', 'Removed nb', 'Removed abundance', 'Abundance of the most abundant removed', 'Detected nb', 'Detected abundance', 'Abundance of the most abundant detected']) + "\n" )
    for sample in sorted(samples_names):
        sample_remove_results = "\t".join(map(str, [sample,
                                                    log_remove_spl[sample]['nb_kept'],
                                                    log_remove_spl[sample]['kept_abundance'],
                                                    log_remove_spl[sample]['nb_removed'],
                                                    log_remove_spl[sample]['removed_abundance'],
                                                    log_remove_spl[sample]['removed_max_abundance']
                                                    ]))
        FH_out.write( sample_remove_results + "\n" )
    FH_out.write( "\n" )
    FH_out.close()

def remove_itsx_biom(preserve, samples, itsx_file, in_fasta_file, in_biom_file, out_removed, out_biom_file, global_report, bySample_report, log_file ):
    """
    @summary: Removes the chimera observation from BIOM.
    @param preserve: [boolean] preserve sequences and remove only not found sequences
    @param samples: [list] samples name list
    @param chimera_files : [list] samples chimera files
    @param in_biom_file: [str] The path to the BIOM file to filter.
    @param out_biom_file: [str] The path to the BIOM after filter.
    
    @param global_report: [dict] This dictionary is update with the global
                          number of removed observations, the global removed
                          abundance, ...
    @param bySample_report: [dict] This dictionary is update for add by sample the
                            number of removed observations, the removed
                            abundance, ...
    @param log_file : [path] Path to general log output file
    """
    FH_log = Logger(log_file)
    FH_log.write("## Removes the observations after ITSx.\n")
    
    # Init bySample_report
    for sample_name in samples:
        bySample_report[sample_name] = {
            'nb_kept': 0,
            'kept_abundance': 0,
            'nb_removed': 0,
            'removed_abundance': 0,
            'removed_max_abundance': 0
        }
    
    in_biom = BiomIO.from_json(in_biom_file)
    kept_clusters_ids = list()
    FH_out_removed = FastaIO( out_removed, "wt" )
    
    if not preserve:
        ## Retrieve IDs to remove
        # Get initial observation names
        all_clusters_ids = list()
        for observation_name in in_biom.get_observations_counts():
            all_clusters_ids.append(str(observation_name[0]))
        
        # Get kept observation names
        record_iter = FastaIO( itsx_file )
        for idx, record in enumerate(record_iter):
            kept_clusters_ids.append(record.id)
        # Difference between initial and kept
        lost_clusters = [x for x in all_clusters_ids if x not in kept_clusters_ids]
        record_iter = FastaIO( in_fasta_file )
        for record in record_iter:
            if record.id in lost_clusters:
                FH_out_removed.write(record)
                
    else:
        # Retrive IDs to remove
        record_iter = FastaIO( itsx_file )
        lost_clusters = list()
        for idx, record in enumerate(record_iter):
            lost_clusters.append(record.id)
        
        FH_out = FastaIO(itsx_file,"wt")
        record_iter = FastaIO( in_fasta_file )
        for record in record_iter:
            if not record.id in lost_clusters:
                FH_out.write(record)
                kept_clusters_ids.append(record.id)
            else:
                FH_out_removed.write(record)

    FH_out_removed.close()

    # Get abundance metrics of lost observations
    for lost_cluster in lost_clusters:
        global_report['nb_removed'] += 1
        global_report['abundance_removed'] += in_biom.get_observation_count(lost_cluster)
        # By sample metrics
        for sample in in_biom.get_samples_by_observation(lost_cluster):
            sample_count = in_biom.get_count(lost_cluster, sample['id'])
            bySample_report[sample['id']]['nb_removed'] += 1
            bySample_report[sample['id']]['removed_abundance'] += sample_count
            bySample_report[sample['id']]['removed_max_abundance'] = max(bySample_report[sample['id']]['removed_max_abundance'], sample_count)
    
    # Get abundance metrics of kept observations
    for kept_cluster in kept_clusters_ids:
        global_report['nb_kept'] += 1
        global_report['abundance_kept'] += in_biom.get_observation_count(kept_cluster)
        # By sample metrics
        for sample in in_biom.get_samples_by_observation(kept_cluster):
            sample_count = in_biom.get_count(kept_cluster, sample['id'])
            bySample_report[sample['id']]['nb_kept'] += 1
            bySample_report[sample['id']]['kept_abundance'] += sample_count
    
    # Write final output BIOM
    out_biom = BiomIO.from_json(in_biom_file)
    remove_observations( lost_clusters, in_biom_file, out_biom_file )

def main_process(args):
    try : 
        tmpFiles = TmpFiles(os.path.split(args.output_fasta)[0])
        summary = tmpFiles.add(args.summary)
        # set number of process
        nb_seq = get_fasta_nb_seq(args.input_fasta)
        
        in_biom = BiomIO.from_json( args.input_biom )
        samples = list()
        for sample_name in in_biom.get_samples_names():
            samples.append(sample_name)
        
        if args.nb_cpus == 1 or nb_seq < 10:
            in_fasta = os.path.abspath(args.input_fasta)
            tmp_dir = tmpFiles.add_dir(os.path.split(args.output_fasta)[1])
            if not args.check_its_only:
                process_ITSx(in_fasta, args.its, tmp_dir, args.output_fasta, args.log_file)
            else:
                process_ITSx(in_fasta, 'no_detections', tmp_dir, args.output_fasta, args.log_file)
        else:
            fasta_ITSx_list = list()
            ITSx_outputs = list()
            logs_ITSx = list()
            tmp_dirs = list()

            split_fasta(args.input_fasta, tmpFiles, max(1,args.nb_cpus-1), fasta_ITSx_list, args.log_file)
            ITSx_outputs = [tmpFiles.add(os.path.basename(current_fasta) + ".fasta") for current_fasta in fasta_ITSx_list]
            logs_ITSx = [tmpFiles.add(os.path.basename(current_fasta) + "_itsx.log") for current_fasta in fasta_ITSx_list]
            tmp_dirs = [ tmpFiles.add_dir(os.path.split(current_fasta)[1]) for current_fasta in fasta_ITSx_list ]
            if not args.check_its_only:
                parallel_submission( process_ITSx, fasta_ITSx_list, args.its, tmp_dirs, ITSx_outputs, logs_ITSx, len(fasta_ITSx_list) )
            else:
                parallel_submission( process_ITSx, fasta_ITSx_list, 'no_detections', tmp_dirs, ITSx_outputs, logs_ITSx, len(fasta_ITSx_list) )

            # Logs
            append_results(ITSx_outputs, logs_ITSx, args.output_fasta, args.log_file)
        
        log_remove_global = { 'nb_kept': 0,
                              'abundance_kept': 0,
                              'nb_removed': 0,
                              'abundance_removed': 0,
                            }
        log_remove_spl = {}
        if args.input_biom is not None:
            remove_itsx_biom(args.check_its_only, samples, args.output_fasta, args.input_fasta, args.input_biom, args.out_removed, args.output_biom, log_remove_global, log_remove_spl, args.log_file )
        
        # Summary
        write_summary( samples, log_remove_global, log_remove_spl, summary )
        
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
    parser.add_argument( '-i', '--its', type=str, required=True, choices=['ITS1','ITS2'], help='Which ITS region are targeted. either ITS1 or ITS2 ')
    parser.add_argument( '--check-its-only', action='store_true', default=False, help='Check only if sequences seem to be an ITS. No sequence trimming will happen' )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-f', '--input-fasta', required=True, help='The fasta input sequences to treat' )
    group_input.add_argument( '-b', '--input-biom', required=True, help='The abundance file for clusters by sample (format: BIOM).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-fasta', default='ITSX.fasta', help='Fasta with ITS sequences only. [Default: %(default)s]' )
    group_output.add_argument( '-m', '--out-removed', default='removed.fasta', help='sequences file removed (format: fasta). [Default: %(default)s]')
    group_output.add_argument( '--summary', default='summary.tsv', help='Summary file. [Default: %(default)s]' )
    group_output.add_argument( '-a', '--output-biom', default='ITSX.biom', help='Updated BIOM file. [Default: %(default)s]' )
    group_output.add_argument( '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.' )
    args = parser.parse_args()

    # Process
    main_process(args)
