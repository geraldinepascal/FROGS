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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse - Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.7.1'
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

from frogsBiom import BiomIO
from frogsSequenceIO import *
from frogsUtils import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def get_sample_resuts( log_file ):
    """
    @summary: Returns the sample results (number of sequences after each filters).
    @param log_file: [str] Path to a log file.
    @return: [list] The number of sequences after each filter.
    """
    res_dic=dict()
    
    FH_input = open(log_file)
    for line in FH_input:
        if line.strip().startswith('nb_chimera: '):
            res_dic["nb_chimera"]= int(line.split(':')[1].strip())
        elif line.strip().startswith('chimera_abun: '):
            res_dic["chimera_abundance"] = int(line.split(':')[1].strip())
        elif line.strip().startswith('max_chimera_abun: '):
            res_dic["chimera_max_abundance"] = int(line.split(':')[1].strip())
    FH_input.close()
    return res_dic
    
def write_summary( samples_names, sample_logs, log_remove_global, log_remove_spl, out_file ):
    """
    @summary: Writes the summary file.
    @param samples_names: [list] The samples names.
    @param sample_logs: [list] list of sample logs files
    @param log_remove_global: [dict] The global remove metrics.
    @param log_remove_spl: [dict] The remove metrics by sample.
    @param out_file: [str] Path to the summary file.
    """
    # Collect metrics
    detection_results = dict()
    for idx, sample in enumerate(samples_names):
        detection_results[sample] = get_sample_resuts(sample_logs[idx])

    # Writes output
    FH_out = open(out_file, "wt")

    global_remove_results = [ log_remove_global['nb_removed'], log_remove_global['nb_kept'],
                              log_remove_global['abundance_removed'], log_remove_global['abundance_kept'],
                              log_remove_global['nb_ambiguous'], log_remove_global['abundance_ambiguous'] ]
    FH_out.write( '##Metrics global\n' )
    FH_out.write( "\t".join(['#Nb removed', 'Nb kept', 'Abundance removed', 'Abundance kept', 'Nb ambiguous', 'Abundance ambiguous']) + "\n" )
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
                                                    log_remove_spl[sample]['removed_max_abundance'],
                                                    detection_results[sample]['nb_chimera'],
                                                    detection_results[sample]['chimera_abundance'],
                                                    detection_results[sample]['chimera_max_abundance'],
                                                    ]))
        FH_out.write( sample_remove_results + "\n" )
    FH_out.write( "\n" )

    FH_out.close()

def get_obs_from_biom( in_biom ):
    """
    @summary: Returns the counts by observation from a BIOM file.
    @param in_biom: Path to the BIOM.
    @return: [dict] Returns the counts by observation.
    """
    observ_dict = dict()
    biom = BiomIO.from_json(in_biom)
    for observation_name in biom.get_observations_names():
        observ_dict[observation_name] = biom.get_observation_count(observation_name)
    del biom
    return observ_dict

def get_obs_from_count( in_count ):
    """
    @summary: Returns the counts by observation from a COUNT file.
    @param in_count: Path to the COUNT.
    @return: [dict] Returns the counts by observation.
    """
    observ_dict = dict()
    in_count_fh = open(in_count)
    header = in_count_fh.readline()
    for line in in_count_fh:
        line_fields = line.strip().split()
        observation_sum = 0
        for count in line_fields[1:]:
            observation_sum += int(count)
        observ_dict[line_fields[0]] = observation_sum
    in_count_fh.close()
    return observ_dict

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
        raise Exception( "\n\n#ERROR : " + error_msg + "\n\n" )

def remove_chimera_fasta( in_fasta, out_fasta, kept_observ, user_size_separator ):
    in_fasta_fh = FastaIO( in_fasta )
    out_fasta_fh = FastaIO( out_fasta, "wt")
    for record in in_fasta_fh:
        real_id = record.id
        if user_size_separator is not None and user_size_separator in record.id:
            real_id = record.id.rsplit(user_size_separator, 1)[0]
        if real_id in kept_observ:
            record.id = real_id
            if user_size_separator is not None:
                record.id = real_id + user_size_separator + str(kept_observ[real_id])
            out_fasta_fh.write(record)
    in_fasta_fh.close()
    out_fasta_fh.close()

def remove_chimera_biom( samples, chimera_files, in_biom_file, out_biom_file, lenient_filter, global_report, bySample_report, log_file ):
    """
    @summary: Removes the chimera observation from BIOM.
    @param samples: [list] samples name list
    @param chimera_files : [list] samples chimera files
    @param in_biom_file: [str] The path to the BIOM file to filter.
    @param out_biom_file: [str] The path to the BIOM after filter.
    @param lenient_filter: [bool] True: removes one sequence in all samples
                           only if it is detected as chimera in all samples
                           where it is present. With False removes one
                           sequence in all samples if it is detected as chimera
                           in at least one sample.
    @param global_report: [dict] This dictionary is update with the global
                          number of removed observations, the global removed
                          abundance, ...
    @param bySample_report: [dict] This dictionary is update for add by sample the
                            number of removed observations, the removed
                            abundance, ...
    @param log_file : [path] Path to general log output file
    """
    FH_log = Logger(log_file)
    FH_log.write("## Removes the chimera observation from BIOM.\n")
    nb_sample_by_chimera = dict()

    # Init bySample_report
    for sample_name in samples:
        bySample_report[sample_name] = {
            'nb_kept': 0,
            'kept_abundance': 0,
            'nb_removed': 0,
            'removed_abundance': 0,
            'removed_max_abundance': 0
        }

    # Retrieve chimera
    for chimera_file in chimera_files:
        chimera_fh = open( chimera_file)
        for line in chimera_fh:
            observation_name = line.strip()
            if observation_name not in nb_sample_by_chimera:
                nb_sample_by_chimera[observation_name] = 0
            nb_sample_by_chimera[observation_name] += 1
        chimera_fh.close()

    # Remove chimera
    removed_chimera = list()
    biom = BiomIO.from_json(in_biom_file)
    for chimera_name in list(nb_sample_by_chimera.keys()):
        is_always_chimera = True
        nb_sample_with_obs = sum( 1 for sample in biom.get_samples_by_observation(chimera_name) )
        observation_abundance = biom.get_observation_count(chimera_name)
        if nb_sample_with_obs != nb_sample_by_chimera[chimera_name]:
            is_always_chimera = False
            global_report['nb_ambiguous'] += 1
            global_report['abundance_ambiguous'] += observation_abundance
            FH_log.write("'" + chimera_name + "' is not interpreted as chimera in all samples where it is present.\n")
        if not lenient_filter or is_always_chimera:
            removed_chimera.append(chimera_name)
            # Global metrics
            global_report['nb_removed'] += 1
            global_report['abundance_removed'] += observation_abundance
            # By sample metrics
            for sample in biom.get_samples_by_observation(chimera_name):
                sample_count = biom.get_count(chimera_name, sample['id'])
                bySample_report[sample['id']]['nb_removed'] += 1
                bySample_report[sample['id']]['removed_abundance'] += sample_count
                bySample_report[sample['id']]['removed_max_abundance'] = max(bySample_report[sample['id']]['removed_max_abundance'], sample_count)
    biom.remove_observations(removed_chimera)

    # Nb non-chimera
    for observation_name in biom.get_observations_names():
        global_report['nb_kept'] += 1
        global_report['abundance_kept'] += biom.get_observation_count(observation_name)
        # By sample metrics
        for sample in biom.get_samples_by_observation(observation_name):
            sample_count = biom.get_count(observation_name, sample['id'])
            bySample_report[sample['id']]['nb_kept'] += 1
            bySample_report[sample['id']]['kept_abundance'] += sample_count
    BiomIO.write(out_biom_file, biom)
    FH_log.close()

def remove_chimera_count( samples, chimera_files, in_count_file, out_count_file, lenient_filter, global_report, bySample_report, log_file ):
    """
    @summary: Removes the chimera observation from TSV.
    @param samples: [list] samples name list
    @param chimera_files : [list] samples chimera files
    @param in_count_file: [str] The path to the COUNT file to filter.
    @param out_count_file: [str] The path to the COUNT after filter.
    @param lenient_filter: [bool] True: removes one sequence in all samples
                           only if it is detected as chimera in all samples
                           where it is present. With False removes one
                           sequence in all samples if it is detected as chimera
                           in at least one sample.
    @param global_report: [dict] This dictionary is update with the global
                          number of removed observations, the global removed
                          abundance, ...
    @param bySample_report: [dict] This dictionary is update for add by sample the
                            number of removed observations, the removed
                            abundance, ...
    @param log_file : [path] Path to general log output file
    """
    FH_log = Logger(log_file)
    FH_log.write("Removes the chimera observation from TSV.\n")
    chimera = dict()
    # Retrieve chimera
    for idx, sample_name in enumerate(samples):
        chimera_fh = open( chimera_files[idx] )
        for line in chimera_fh:
            observation_name = line.strip()
            if observation_name not in chimera:
                chimera[observation_name] = dict()
            chimera[observation_name][sample_name] = True
        chimera_fh.close()

    # Remove chimera
    in_count_fh = open( in_count_file )
    out_count_fh = open( out_count_file, "wt" )
    samples_pos = dict()
    #    header
    header = in_count_fh.readline()
    out_count_fh.write(header)
    for idx, sample_name in enumerate(header.strip().split()[1:]):
        samples_pos[sample_name] = idx
        if sample_name not in bySample_report:
            bySample_report[sample_name] = {
                'nb_kept': 0,
                'kept_abundance': 0,
                'nb_removed': 0,
                'removed_abundance': 0,
                'removed_max_abundance': 0
            }
    #    body
    for line in in_count_fh:
        line_fields = line.strip().split()
        observation_name = line_fields[0]
        observation_counts = [int(sample_count) for sample_count in line_fields[1:]]
        if observation_name not in chimera:
            out_count_fh.write( line )
            global_report['nb_kept'] += 1
            global_report['abundance_kept'] += sum(observation_counts)
            # By sample metrics
            for sample_name in list(samples_pos.keys()):
                sample_count = int(observation_counts[samples_pos[sample_name]])
                if sample_count > 0:
                    bySample_report[sample_name]['nb_kept'] += 1
                    bySample_report[sample_name]['kept_abundance'] += sample_count
        else: # is chimera in at least one sample
            is_always_chimera = True
            for sample_name in list(samples_pos.keys()):
                if sample_name not in chimera[observation_name] and int(observation_counts[samples_pos[sample_name]]) != 0:
                    is_always_chimera = False
            if not is_always_chimera: # is not chimera in all samples where it is find
                global_report['nb_ambiguous'] += 1
                global_report['abundance_ambiguous'] += sum(observation_counts)
                FH_log.write( "'" + observation_name + "' is not interpreted as chimera in all samples where it is present.\n")
            if is_always_chimera or not lenient_filter:
                global_report['nb_removed'] += 1
                global_report['abundance_removed'] += sum(observation_counts)
                # By sample metrics
                for sample_name in list(samples_pos.keys()):
                    sample_count = int(observation_counts[samples_pos[sample_name]])
                    if sample_count > 0:
                        bySample_report[sample_name]['nb_removed'] += 1
                        bySample_report[sample_name]['removed_abundance'] += sample_count
                        bySample_report[sample_name]['removed_max_abundance'] = max(bySample_report[sample_name]['removed_max_abundance'], sample_count)
            else:
                global_report['nb_kept'] += 1
                global_report['abundance_kept'] += sum(observation_counts)
                out_count_fh.write( line )
                # By sample metrics
                for sample_name in list(samples_pos.keys()):
                    sample_count = int(observation_counts[samples_pos[sample_name]])
                    if sample_count > 0:
                        bySample_report[sample_name]['nb_kept'] += 1
                        bySample_report[sample_name]['kept_abundance'] += sample_count

    in_count_fh.close()
    out_count_fh.close()
    FH_log.close()

def chimera( sample_names, input_fasta, input_abund, outputs_fasta, outputs_chimera, log_chimera, user_size_separator ):
    for idx in range(len(sample_names)):
        chimera_by_sample( sample_names[idx], input_fasta, input_abund, outputs_fasta[idx], outputs_chimera[idx], log_chimera[idx], user_size_separator )
    time.sleep(0.5) # Wait to fix 'Exception in thread QueueFeederThread' in slow systems 

def chimera_by_sample( sample_name, input_fasta, input_abund, output_fasta, output_chimera, log_chimera, user_size_separator ):
    tmp_fasta = output_fasta + ".tmp"
    tmp_log = output_fasta + ".log"
    tmp_stderr = output_fasta + ".stderr"
    tmp_stdout = output_fasta + ".stdout"
    size_separator = ";size="
    count_by_obs = dict()

    try:
        FH_log = open(log_chimera,"wt")
        FH_log.write("##Sample : " + sample_name + "\n")
        # Get count by obs
        in_obs_fh = open( input_abund )
        header_line = in_obs_fh.readline().strip()
        sample_idx = header_line.split().index( sample_name )
        for line in in_obs_fh:
            line_fields = line.strip().split()
            if int(line_fields[sample_idx]) != 0:
                count_by_obs[line_fields[0]] = line_fields[sample_idx]
        in_obs_fh.close()

        # Write fasta with observation size
        nb_seq_sample = 0
        in_fasta_fh = FastaIO( input_fasta )
        tmp_fasta_fh = FastaIO( tmp_fasta, "wt")
        for record in in_fasta_fh:
            real_id = record.id
            if user_size_separator is not None and user_size_separator in record.id:
                real_id = record.id.rsplit(user_size_separator, 1)[0]
            if real_id in count_by_obs:
                nb_seq_sample += 1
                record.id = real_id + size_separator + str(count_by_obs[real_id])
                tmp_fasta_fh.write(record)
        in_fasta_fh.close()
        tmp_fasta_fh.close()

        # Chimera cleanning
        if nb_seq_sample != 0:
            FH_log.write("## Vsearch command: " + " ".join(["vsearch", "--uchime_denovo", tmp_fasta, "--nonchimeras", output_fasta, "--uchimeout", tmp_log]) + "\n" )
            submit_cmd( ["vsearch", "--uchime_denovo", tmp_fasta, "--nonchimeras", output_fasta, "--uchimeout", tmp_log], tmp_stdout, tmp_stderr )
        else: # The sample is empty
            FH_log.write("## Empty sample, no chimera research\n")
            open( output_fasta, "wt" ).close()
            open( tmp_log, "wt" ).close()

        # Log
        nb_chimera = 0
        chimera_abun = 0
        max_chimera_abun = 0
        nb_non_chimera = 0
        non_chimera_abun = 0
        in_log_fh = open( tmp_log )
        out_chimera_fh = open( output_chimera, "wt" )
        for line in in_log_fh:
            line_fields = line.strip().split()
            observation_name, size = line_fields[1].rsplit( size_separator, 1 )
            size = int(size)
            if line_fields[-1] == "Y":
                out_chimera_fh.write( observation_name + "\n" )
                nb_chimera += 1
                chimera_abun += size
                max_chimera_abun = max(max_chimera_abun, size)
            else:
                nb_non_chimera += 1
                non_chimera_abun += size
        in_log_fh.close()
        out_chimera_fh.close()
        
        FH_log.write("##Results\n")
        FH_log.write("sample_name: " + sample_name + "\n" + \
                     "nb_chimera: " + str(nb_chimera)  + "\n" + \
                     "chimera_abun: " + str(chimera_abun)  + "\n" + \
                     "max_chimera_abun: " + str(max_chimera_abun)  + "\n" + \
                     "nb_non_chimera: " + str(nb_non_chimera) + "\n" + \
                     "non_chimera_abun: " + str(non_chimera_abun) + "\n" )
        FH_log.close()       
    finally:
        if os.path.exists(tmp_fasta): os.remove(tmp_fasta)
        if os.path.exists(tmp_log): os.remove(tmp_log)
        if os.path.exists(tmp_stdout): os.remove(tmp_stdout)
        if os.path.exists(tmp_stderr): os.remove(tmp_stderr)
        
def log_append_files( log_file, appended_files ):
    """
    @summary: Append content of several log files in one log file.
    @param log_file: [str] The log file where contents of others are appended.
    @param appended_files: [list] List of log files to append.
    """
    FH_log = Logger( log_file )
    FH_log.write( "\n" )
    for current_file in appended_files:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write( line )
        FH_input.close()
        FH_log.write( "\n" )
    FH_log.write( "\n" )
    FH_log.close()
    
def main_process(args):
    tmp_files = TmpFiles(os.path.split(args.non_chimera)[0])

    try:
        if args.out_abundance is None:
            args.out_abundance = "count.tsv"
            if args.biom is not None:
                args.out_abundance = "abundance.biom"

        count_table = args.count
        if args.biom is not None:
            count_table = tmp_files.add("tmp_count.tsv")
            biom = BiomIO.from_json( args.biom )
            BiomIO.write_count_table( count_table, biom )
            del biom

        # Get samples
        samples = list()
        fasta_files=list()
        chimera_files=list()
        sample_logs=list()
        in_count_fh = open( count_table )
        header_line = in_count_fh.readline().strip()
        for sample_name in header_line.split()[1:]:
            samples.append(sample_name)
            fasta_files.append(tmp_files.add(sample_name + ".fasta"))
            chimera_files.append(tmp_files.add(sample_name + ".chimera"))
            sample_logs.append(tmp_files.add(sample_name + "_log.txt"))
        in_count_fh.close()

        # Find chimera
        nb_processses_used = min( len(samples), args.nb_cpus )
        processes = [{'process':None, 'in_file':[], 'out_file':[], 'sample_name':[], 'log' :[]} for idx in range(nb_processses_used)]
        #    Set processes
        for idx, sample_name in enumerate(samples):
            process_idx = idx % nb_processses_used
            processes[process_idx]['sample_name'].append( sample_name )
            processes[process_idx]['in_file'].append( fasta_files[idx] )
            processes[process_idx]['out_file'].append( chimera_files[idx] )
            processes[process_idx]['log'].append( sample_logs[idx] )
        #    Launch processes
        for current_process in processes:
            if idx == 0: # First process is threaded with parent job
                current_process['process'] = threading.Thread( target=chimera, 
                                                               args=(current_process['sample_name'], args.sequences, count_table, current_process['in_file'], current_process['out_file'], current_process['log'], args.size_separator) )
            else: # Others processes are processed on different CPU
                current_process['process'] = multiprocessing.Process( target=chimera, 
                                                                      args=(current_process['sample_name'], args.sequences, count_table, current_process['in_file'], current_process['out_file'], current_process['log'], args.size_separator) )
            current_process['process'].start()
        #    Wait processes end
        for current_process in processes:
            current_process['process'].join()
        #    Check processes status
        for current_process in processes:
            if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
                sys.exit(1)

        # Append independant log files
        log_append_files( args.log_file, sample_logs )

        # Remove chimera
        log_remove_global = { 'nb_kept': 0,
                              'abundance_kept': 0,
                              'nb_removed': 0,
                              'abundance_removed': 0,
                              'nb_ambiguous': 0,
                              'abundance_ambiguous': 0}
        log_remove_spl = {}

        if args.biom is not None:
            remove_chimera_biom( samples, chimera_files, args.biom, args.out_abundance, args.lenient_filter, log_remove_global, log_remove_spl, args.log_file )
            remove_chimera_fasta( args.sequences, args.non_chimera, get_obs_from_biom(args.out_abundance), args.size_separator )
        else:
            remove_chimera_count( samples, chimera_files, args.count, args.out_abundance, args.lenient_filter, log_remove_global, log_remove_spl, args.log_file )
            remove_chimera_fasta( args.sequences, args.non_chimera, get_obs_from_count(args.out_abundance), args.size_separator )

        # Summary
        write_summary( samples, sample_logs, log_remove_global, log_remove_spl, args.summary )
        
    finally:
        if not args.debug:
            tmp_files.deleteAll()

def get_vsearch_version():
    """
    @summary: Return the vserach version.
    @return: [str] The vsearch version.
    """
    version = None
    try:
        cmd = ["vsearch", "--version"]
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        version = stderr.split(",")[0].split()[1] # vsearch v1.1.3_linux_x86_64, 126.0GB RAM, 32 cores           
    except:
        version = "unknown"
    return version


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Removes PCR chimera by samples.'
    )
    parser.add_argument( '--size-separator', help="The size separator if the cluster IDs contain the number of represented sequence (format: '<ID_IN_ABUND_FILE><size_separator><NB_SEQ>'" )
    parser.add_argument( '-l', '--lenient-filter', default=False, action='store_true', help="Removes one sequence in all samples only if it is detected as chimera in all samples where it is present. Without this option the program removes one sequence in all samples if it is detected as chimera in at least one sample." )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ + " [vsearch " + get_vsearch_version() + "]" )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-s', '--sequences', required=True, help='The cluster sequences.' )
    group_exclusion_abundance = group_input.add_mutually_exclusive_group()
    group_exclusion_abundance.add_argument( '-b', '--biom', help='The abundance file for clusters by sample (format: BIOM).' )
    group_exclusion_abundance.add_argument( '-c', '--count', help='The abundance file for clusters by sample (format: count).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-n', '--non-chimera', default='non_chimera.fasta', help='Fasta without chimera. [Default: %(default)s]' )
    group_output.add_argument( '-a', '--out-abundance', default=None, help='Abundance file without chimera.' )
    group_output.add_argument( '--summary', default='summary.tsv', help='Summary file. [Default: %(default)s]' )
    group_output.add_argument( '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.' )
    args = parser.parse_args()

    # Process
    main_process(args)
