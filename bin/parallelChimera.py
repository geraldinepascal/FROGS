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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.5.1'
__email__ = 'frogs@toulouse.inra.fr'
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
from biom import *
from sequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
class TmpFiles:
    """
    @summary: Manager for temporary files.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)

def write_summary( samples_names, log_detection, log_remove_global, log_remove_spl, out_file ):
    """
    @summary: Writes the summary file.
    @param samples_names: [list] The samples names.
    @param log_detection: [Queue] The chimera detection metrics by individual sample.
    @param log_remove_global: [dict] The global remove metrics.
    @param log_remove_spl: [dict] The remove metrics by sample.
    @param out_file: [str] Path to the summary file.
    """
    # Collect metrics
    detection_results = dict()
    for elt in samples_names:
        sample_info = log_detection.get()
        detection_results[sample_info['sample_name']] = {
            'nb_chimera': sample_info['nb_chimera'],
            'chimera_abundance': sample_info['chimera_abun'],
            'chimera_max_abundance': sample_info['max_chimera_abun']
        }

    # Writes output
    FH_out = open(out_file, "w")

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

def remove_chimera_fasta( in_fasta, out_fasta, kept_observ, user_size_separator ):
    in_fasta_fh = FastaIO( in_fasta )
    out_fasta_fh = FastaIO( out_fasta, "w")
    for record in in_fasta_fh:
        real_id = record.id
        if user_size_separator is not None and user_size_separator in record.id:
            real_id = record.id.rsplit(user_size_separator, 1)[0]
        if kept_observ.has_key(real_id):
            record.id = real_id
            if user_size_separator is not None:
                record.id = real_id + user_size_separator + str(kept_observ[real_id])
            out_fasta_fh.write(record)
    in_fasta_fh.close()
    out_fasta_fh.close()

def remove_chimera_biom( samples, in_biom_file, out_biom_file, lenient_filter, global_report, bySample_report ):
    """
    @summary: Removes the chimera observation from BIOM.
    @param samples: [dict] The chimera observations by sample. Example for
                    sample splA: sample['splA']['chimera_path'] where the value
                    is the path to the file containing the list of the chimera
                    observations names.
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
    """
    nb_sample_by_chimera = dict()

    # Init bySample_report
    for sample_name in samples.keys():
        bySample_report[sample_name] = {
            'nb_kept': 0,
            'kept_abundance': 0,
            'nb_removed': 0,
            'removed_abundance': 0,
            'removed_max_abundance': 0
        }

    # Retrieve chimera
    for sample_name in samples.keys():
        chimera_fh = open( samples[sample_name]['chimera_path'] )
        for line in chimera_fh:
            observation_name = line.strip()
            if not nb_sample_by_chimera.has_key(observation_name):
                nb_sample_by_chimera[observation_name] = 0
            nb_sample_by_chimera[observation_name] += 1
        chimera_fh.close()

    # Remove chimera
    removed_chimera = list()
    biom = BiomIO.from_json(in_biom_file)
    for chimera_name in nb_sample_by_chimera.keys():
        is_always_chimera = True
        nb_sample_with_obs = sum( 1 for sample in biom.get_samples_by_observation(chimera_name) )
        observation_abundance = biom.get_observation_count(chimera_name)
        if nb_sample_with_obs != nb_sample_by_chimera[chimera_name]:
            is_always_chimera = False
            global_report['nb_ambiguous'] += 1
            global_report['abundance_ambiguous'] += observation_abundance
            print "'" + chimera_name + "' is not interpreted as chimera in all samples where it is present."
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

def remove_chimera_count( samples, in_count_file, out_count_file, lenient_filter, global_report, bySample_report ):
    """
    @summary: Removes the chimera observation from TSV.
    @param samples: [dict] The chimera observations by sample. Example for
                    sample splA: sample['splA']['chimera_path'] where the value
                    is the path to the file containing the list of the chimera
                    observations names.
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
    """
    chimera = dict()

    # Retrieve chimera
    for sample_name in samples.keys():
        chimera_fh = open( samples[sample_name]['chimera_path'] )
        for line in chimera_fh:
            observation_name = line.strip()
            if not chimera.has_key(observation_name):
                chimera[observation_name] = dict()
            chimera[observation_name][sample_name] = True
        chimera_fh.close()

    # Remove chimera
    in_count_fh = open( in_count_file )
    out_count_fh = open( out_count_file, "w" )
    samples_pos = dict()
    #    header
    header = in_count_fh.readline()
    out_count_fh.write(header)
    for idx, sample_name in enumerate(header.strip().split()[1:]):
        samples_pos[sample_name] = idx
        if not bySample_report.has_key( sample_name ):
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
        if not chimera.has_key(observation_name):
            out_count_fh.write( line )
            global_report['nb_kept'] += 1
            global_report['abundance_kept'] += sum(observation_counts)
            # By sample metrics
            for sample_name in samples_pos.keys():
                sample_count = int(observation_counts[samples_pos[sample_name]])
                if sample_count > 0:
                    bySample_report[sample_name]['nb_kept'] += 1
                    bySample_report[sample_name]['kept_abundance'] += sample_count
        else: # is chimera in at least one sample
            is_always_chimera = True
            for sample_name in samples_pos.keys():
                if not chimera[observation_name].has_key(sample_name) and int(observation_counts[samples_pos[sample_name]]) != 0:
                    is_always_chimera = False
            if not is_always_chimera: # is not chimera in all samples where it is find
                global_report['nb_ambiguous'] += 1
                global_report['abundance_ambiguous'] += sum(observation_counts)
                print "'" + observation_name + "' is not interpreted as chimera in all samples where it is present."
            if is_always_chimera or not lenient_filter:
                global_report['nb_removed'] += 1
                global_report['abundance_removed'] += sum(observation_counts)
                # By sample metrics
                for sample_name in samples_pos.keys():
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
                for sample_name in samples_pos.keys():
                    sample_count = int(observation_counts[samples_pos[sample_name]])
                    if sample_count > 0:
                        bySample_report[sample_name]['nb_kept'] += 1
                        bySample_report[sample_name]['kept_abundance'] += sample_count

    in_count_fh.close()
    out_count_fh.close()

def chimera( sample_names, input_fasta, input_abund, outputs_fasta, outputs_chimera, log_chimera, user_size_separator ):
    for idx in range(len(sample_names)):
        chimera_by_sample( sample_names[idx], input_fasta, input_abund, outputs_fasta[idx], outputs_chimera[idx], log_chimera, user_size_separator )
    time.sleep(0.5) # Wait to fix 'Exception in thread QueueFeederThread' in slow systems 

def chimera_by_sample( sample_name, input_fasta, input_abund, output_fasta, output_chimera, log_chimera, user_size_separator ):
    tmp_fasta = output_fasta + ".tmp"
    tmp_log = output_fasta + ".log"
    tmp_stderr = output_fasta + ".stderr"
    tmp_stdout = output_fasta + ".stdout"
    size_separator = ";size="
    count_by_obs = dict()

    try:
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
        tmp_fasta_fh = FastaIO( tmp_fasta, "w")
        for record in in_fasta_fh:
            real_id = record.id
            if user_size_separator is not None and user_size_separator in record.id:
                real_id = record.id.rsplit(user_size_separator, 1)[0]
            if count_by_obs.has_key(real_id):
                nb_seq_sample += 1
                record.id = real_id + size_separator + str(count_by_obs[real_id])
                tmp_fasta_fh.write(record)
        in_fasta_fh.close()
        tmp_fasta_fh.close()

        # Chimera cleanning
        if nb_seq_sample != 0:
            submit_cmd( ["vsearch", "--uchime_denovo", tmp_fasta, "--nonchimeras", output_fasta, "--uchimeout", tmp_log], tmp_stdout, tmp_stderr )
        else: # The sample is empty
            open( output_fasta, "w" ).close()
            open( tmp_log, "w" ).close()

        # Log
        nb_chimera = 0
        chimera_abun = 0
        max_chimera_abun = 0
        nb_non_chimera = 0
        non_chimera_abun = 0
        in_log_fh = open( tmp_log )
        out_chimera_fh = open( output_chimera, "w" )
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
        log_chimera.put({ 'sample_name': sample_name,
                          'nb_chimera': nb_chimera,
                          'chimera_abun': chimera_abun,
                          'max_chimera_abun': max_chimera_abun,
                          'nb_non_chimera': nb_non_chimera,
                          'non_chimera_abun': non_chimera_abun
        })
    finally:
        if os.path.exists(tmp_fasta): os.remove(tmp_fasta)
        if os.path.exists(tmp_log): os.remove(tmp_log)
        if os.path.exists(tmp_stdout): os.remove(tmp_stdout)
        if os.path.exists(tmp_stderr): os.remove(tmp_stderr)

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
        samples = dict()
        in_count_fh = open( count_table )
        header_line = in_count_fh.readline().strip()
        for sample_name in header_line.split()[1:]:
            samples[sample_name] = { 'fasta_path': tmp_files.add(sample_name + ".fasta"),
                                     'chimera_path': tmp_files.add(sample_name + ".chimera")
            }
        in_count_fh.close()

        # Find chimera
        log_detection = Queue()
        nb_processses_used = min( len(samples.keys()), args.nb_cpus )
        processes = [{'process':None, 'in_file':[], 'out_file':[], 'sample_name':[]} for idx in range(nb_processses_used)]
        #    Set processes
        for idx, sample_name in enumerate(samples.keys()):
            process_idx = idx % nb_processses_used
            processes[process_idx]['sample_name'].append( sample_name )
            processes[process_idx]['in_file'].append( samples[sample_name]['fasta_path'] )
            processes[process_idx]['out_file'].append( samples[sample_name]['chimera_path'] )
        #    Launch processes
        for current_process in processes:
            if idx == 0: # First process is threaded with parent job
                current_process['process'] = threading.Thread( target=chimera, 
                                                               args=(current_process['sample_name'], args.sequences, count_table, current_process['in_file'], current_process['out_file'], log_detection, args.size_separator) )
            else: # Others processes are processed on different CPU
                current_process['process'] = multiprocessing.Process( target=chimera, 
                                                                      args=(current_process['sample_name'], args.sequences, count_table, current_process['in_file'], current_process['out_file'], log_detection, args.size_separator) )
            current_process['process'].start()
        #    Wait processes end
        for current_process in processes:
            current_process['process'].join()
        #    Check processes status
        for current_process in processes:
            if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
                sys.exit(1)

        # Remove chimera
        log_remove_global = { 'nb_kept': 0,
                              'abundance_kept': 0,
                              'nb_removed': 0,
                              'abundance_removed': 0,
                              'nb_ambiguous': 0,
                              'abundance_ambiguous': 0}
        log_remove_spl = {}

        if args.biom is not None:
            remove_chimera_biom( samples, args.biom, args.out_abundance, args.lenient_filter, log_remove_global, log_remove_spl )
            remove_chimera_fasta( args.sequences, args.non_chimera, get_obs_from_biom(args.out_abundance), args.size_separator )
        else:
            remove_chimera_count( samples, args.count, args.out_abundance, args.lenient_filter, log_remove_global, log_remove_spl )
            remove_chimera_fasta( args.sequences, args.non_chimera, get_obs_from_count(args.out_abundance), args.size_separator )

        # Summary
        write_summary( samples.keys(), log_detection, log_remove_global, log_remove_spl, args.summary )
    finally:
        if not args.debug:
            tmp_files.deleteAll()


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
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used." )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-s', '--sequences', required=True, help='The cluster sequences.' )
    group_exclusion_abundance = group_input.add_mutually_exclusive_group()
    group_exclusion_abundance.add_argument( '-b', '--biom', help='The abundance file for clusters by sample (format: BIOM).' )
    group_exclusion_abundance.add_argument( '-c', '--count', help='The abundance file for clusters by sample (format: count).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-n', '--non-chimera', default='non_chimera.fasta', help='Fasta without chimera.')
    group_output.add_argument( '-a', '--out-abundance', default=None, help='Abundance file without chimera.')
    group_output.add_argument( '--summary', default='summary.tsv', help='Summary file.')
    args = parser.parse_args()

    # Process
    main_process(args)
