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

__author__ = 'Katia Vidal - Team NED Toulouse AND Frederic Escudie - Plateforme bioinformatique Toulouse AND Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '4.1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'


import os
import sys
import json
import operator
import argparse
import collections

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

from frogsUtils import *
from frogsBiom import Biom, BiomIO
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class Depths(Cmd):
    """
    @summary: Writes by abundance the number of clusters.
    """
    def __init__(self, in_biom, out_tsv):
        """
        @param in_biom: [str] The processed BIOM path.
        @param out_tsv: [str] The path of the output.
        """
        Cmd.__init__( self,
                      'biomTools.py',
                      'Writes by abundance the number of clusters.',
                      'obsdepth --input-file ' + in_biom + ' --output-file ' + out_tsv,
                      '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip() 

class UpdateFasta(Cmd):
    """
    @summary: Updates fasta file based on sequence in biom file
    """
    def __init__(self, in_biom, in_fasta, out_fasta, log):
        """
        @param in_biom: [str] Path to BIOM file.
        @param nb_read : [int] Number of reads per sample
        @param out_biom: [str] Path to output BIOM file.
        """
        Cmd.__init__( self,
                      'biomFastaUpdate.py',
                      'Updates fasta file based on sequence in biom file.',
                      "--input-biom " + in_biom + " --input-fasta " + in_fasta + " --output-file " + out_fasta + " --log " + log,
                      '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()                      

class RemoveConta(Cmd):
    def __init__(self, contaminant_db, in_fasta, in_biom, cleaned_fasta, cleaned_biom, contaminated_fasta, log, nb_cpus, debug ):
        """
        @param contaminant_db: [str] Path to the databank.
        @param in_fasta: [str] Path to the sequences file.
        @param in_biom: [str] Path to the abundance file.
        @param cleaned_fasta: [str] Path to the cleaned sequences file.
        @param cleaned_biom: [str] Path to the cleaned abundance file.
        @param contaminated_fasta: [str] Path to the contaminated sequences file.
        @param log: [str] Path to the log file.
        @param nb_cpus: [int] Number of usable CPUs.
        """
        command_opt = "--nb-cpus " + str(nb_cpus) + " --word-size 40 --min-identity 0.8 --min-coverage 0.8 --input-fasta " + in_fasta + " --contaminant-db " + contaminant_db + " --input-biom " + in_biom + " --clean-fasta " + cleaned_fasta + " --clean-biom " + cleaned_biom + " --conta-fasta " + contaminated_fasta + " --log-file " + log
        if debug:
            command_opt += " --debug"
        Cmd.__init__( self,
                      "removeConta.py",
                      "Removes contaminant sequences.",
                      command_opt,
                      "--version")
        self.process_log = log
        self.contaminated_fasta = contaminated_fasta

    def write_conta_list(self, excluded_list):
        """
        @summary: Writes in output file the list of the sequences IDs tagged as contaminant.
        @excluded_list: [str] The path to the output file.
        """
        FH_contaminated_fasta = FastaIO( self.contaminated_fasta )
        FH_contaminated_list = open( excluded_list, "wt" )
        for record in FH_contaminated_fasta:
            FH_contaminated_list.write( record.id + "\n" )
        FH_contaminated_fasta.close()
        FH_contaminated_list.close()

    def parser(self, log_file):
        """
        @summary: Parse the command results to add information in log_file.
        @log_file: [str] Path to the process log file.
        """
        FH_process_log = open(self.process_log, "rt")
        for line in FH_process_log:
            if line.startswith("#Processed"):
                nb_processed = line.split(":")[1].strip()
            elif line.startswith("#Contaminated"):
                nb_contaminated = line.split(":")[1].strip()
        FH_process_log.close()
        FH_log = Logger(log_file)
        FH_log.write('Results:\n')
        FH_log.write('\tnumber of processed sequences: ' + str(nb_processed) + '\n')
        FH_log.write('\tnumber of removed contaminated sequences: ' + str(nb_contaminated) + '\n\n')
        FH_log.close()
        
    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()        


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def excluded_obs_on_samplePresence(input_biom, min_sample_presence, excluded_file):
    """
    @summary: Writes the list of the observations present in an insufficient number of samples.
    @param input_biom: [str] The path to the BIOM file to check.
    @param min_sample_presence: [int] The observations present in a number of samples inferior than this value are reported in the excluded file.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "wt" )
    for observation_name in biom.get_observations_names():
        nb_samples = sum(1 for x in biom.get_samples_by_observation(observation_name))
        if nb_samples < min_sample_presence:
            FH_excluded_file.write( observation_name + "\n" )
    FH_excluded_file.close()

def excluded_obs_on_replicatePresence(input_biom, replicate_file, min_replicate_presence, log_file, excluded_file):
    """
    @summary:
    """
    # Files handling
    FH_log = Logger(log_file)
    FH_log.write('#Replicate groups:\n')
    biom = BiomIO.from_json( input_biom )
    FH_replicate_file = open(replicate_file).readlines()
    lines = [line for line in FH_replicate_file if line.strip()]
    FH_excluded_file = open( excluded_file, "wt" )
    # Indentify replicates
    groups_to_replicates = collections.defaultdict(list)
    samples_to_search = list()
    for l in lines:
        l = l.strip().split()
        if l[0] not in biom.get_samples_names():
            raise_exception( Exception( "\n\n#ERROR : "+l[0]+" sample not in biom file.\n\n" ))
        groups_to_replicates["".join(l[1:]).strip()].append(l[0])
        samples_to_search.append(l[0])
    # Writes replicate groups into log file
    for group, replicates in groups_to_replicates.items():
        FH_log.write(group+"\t"+", ".join(replicates)+"\n")
    # Search the excluded clusters
    for observation_name in biom.get_observations_names():
        groups_to_counts = collections.defaultdict(int)
        for sample in biom.get_samples_by_observation(observation_name):
            if sample['id'] in samples_to_search:
                group = [group for group, replicates in groups_to_replicates.items() if sample['id'] in replicates][0]
                groups_to_counts[group] += 1
        to_exclude = True
        for group, count in groups_to_counts.items():
            if min_replicate_presence * len(groups_to_replicates[group]) <= count:
                to_exclude = False
        if to_exclude:  
            FH_excluded_file.write( observation_name + "\n" )
    FH_log.close()
    FH_excluded_file.close()

def excluded_obs_on_abundance(input_biom, min_abundance, excluded_file):
    """
    @summary: Writes the list of the observations with an insufficient abundance.
    @param input_biom: [str] The path to the BIOM file to check.
    @param min_abundance: [int/float] The observations with an abundance inferior than this value are reported in the excluded file.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "wt" )
    min_nb_seq = min_abundance
    if type(min_abundance) == float:
        min_nb_seq = biom.get_total_count() * min_abundance
    for idx, count_by_sample in enumerate(biom.to_count()):
        observation = biom.rows[idx]
        abundance = sum(count_by_sample)
        if abundance < min_nb_seq:
            FH_excluded_file.write( str(observation["id"]) + "\n" )
    FH_excluded_file.close()


def excluded_obs_on_nBiggest( input_biom, nb_selected, excluded_file ):
    """
    @summary: Writes the list of all the observations without the n most abundant.
    @param input_biom: [str] The path to the BIOM file.
    @param threshold: [float] The number of the most abundant observations that will not be written in the excluded list.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "wt" )
    sorted_obs_counts = sorted( biom.get_observations_counts(), key=lambda observation: observation[1], reverse=True )
    for observation_name, observation_count in sorted_obs_counts[nb_selected:]:
        FH_excluded_file.write( observation_name + "\n" )
    FH_excluded_file.close()

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

def uniq_from_files_lists( in_files ):
    """
    @summary: Returns an list without duplicated elements from several list files.
    @param in_files: [list] The list of files paths. Each file contains a list.
    @return: [list] The list without duplicated elements.
    """
    uniq = dict()
    for current_file in in_files:
        FH_current_file = open( current_file )
        for line in FH_current_file:
            uniq[line.strip()] = 1
        FH_current_file.close()
    return list(uniq.keys())

def write_exclusion( discards, excluded_file ):
    """
    @summary: Writes the list of observations removed by each filter.
    @param discards: [dict] By filter the path of the file that contains the list of the removed observations.
    @param excluded_file: [str] The path to the output file.
    """
    FH_excluded = open( excluded_file, "wt" )
    list_FH_discards = list()

    # Header
    header_line_fields = list()
    for filter in discards:
        header_line_fields.append( filter )
        list_FH_discards.append( open(discards[filter]) )
    FH_excluded.write( "#" + "\t".join(header_line_fields)  + "\n" )

    # Excluded
    nb_eof = 0
    while nb_eof < len(list_FH_discards):
        discards_line_fields = list()
        for FH_idx, FH_curent_filter in enumerate(list_FH_discards): # For each filter
            observation = ""
            if FH_curent_filter is not None: # Process next line if the discard file is not closed
                observation = FH_curent_filter.readline()
                if observation == "":
                    FH_curent_filter.close()
                    list_FH_discards[FH_idx] = None
                    nb_eof += 1
            discards_line_fields.append( observation.strip() )
        if nb_eof < len(list_FH_discards):
            FH_excluded.write( "\t".join(discards_line_fields)  + "\n" )
    FH_excluded.close()

def write_summary( summary_file, input_biom, output_biom, replicate_log, discards, depth_file ):
    """
    @summary: Writes the process summary.
    @param summary_file: [str] The path to the output file.
    @param input_biom: [str] The path to the BIOM before program execution.
    @param output_biom: [str] The path to the BIOM after program execution.
    @param discards: [dict] By filter the path of the file that contains the list of the removed observations.
    """
    global_results = {
        'nb_clstr_kept': 0,
        'nb_clstr_ini': 0,
        'nb_seq_kept': 0,
        'nb_seq_ini': 0
    }
    samples_results = dict()
    filters_results = dict()

    # Global before filters
    in_biom = BiomIO.from_json( input_biom )
    for observation_name in in_biom.get_observations_names():
        global_results['nb_clstr_ini'] += 1
        global_results['nb_seq_ini'] += in_biom.get_observation_count( observation_name )
    for sample_name in in_biom.get_samples_names():
        samples_results[sample_name] = {
            'initial': sum( 1 for x in in_biom.get_observations_by_sample(sample_name) ),
            'filtered': dict(),
            'kept': 0
        }

    # By sample and by filters
    filters_intersections = dict()
    for filter in list(discards.keys()):
        FH_filter = open( discards[filter] )
        for line in FH_filter:
            observation_name = line.strip()
            if observation_name not in filters_intersections:
                filters_intersections[observation_name] = dict()
            filters_intersections[observation_name][filter] = 1
        FH_filter.close()
    for observation_name in list(filters_intersections.keys()):
        # Removed intersection
        intersections_key = "--@@--".join(sorted( filters_intersections[observation_name].keys() ))
        if intersections_key not in filters_results:
            filters_results[intersections_key] = {
                'filters': list(filters_intersections[observation_name].keys()),
                'count': 0
            }
        filters_results[intersections_key]['count'] += 1

        # Filters by samples
        for sample in in_biom.get_samples_by_observation(observation_name):
            for filter in filters_intersections[observation_name]:
                if filter not in samples_results[sample['id']]['filtered']:
                    samples_results[sample['id']]['filtered'][filter] = 0
                samples_results[sample['id']]['filtered'][filter] += 1
    del in_biom
    
    # Get size distribution data
    clusters_size = list()
    counts = list()
    FH_depth = open( depth_file )
    for line in FH_depth:
        if not line.startswith('#'):
            fields = line.strip().split()
            if fields[1] != "0":
                clusters_size.append( int(fields[0]) )
                counts.append( int(fields[1]) )
    FH_depth.close()

    # Get sample data
    biom = BiomIO.from_json( input_biom )
    samples_distrib = dict()
    for sample_name in biom.get_samples_names():
        shared_seq = 0
        shared_observations = 0
        own_seq = 0
        own_observations = 0
        for observation in biom.get_observations_by_sample(sample_name):
            obs_count_in_spl = biom.get_count( observation['id'], sample_name )
            if obs_count_in_spl != 0 and obs_count_in_spl == biom.get_observation_count(observation['id']):
                own_observations += 1
                own_seq += obs_count_in_spl
            else:
                shared_observations += 1
                shared_seq += obs_count_in_spl
        samples_distrib[sample_name] = {
            'shared_seq': shared_seq,
            'shared_observations': shared_observations,
            'own_seq': own_seq,
            'own_observations': own_observations
        }
    del biom

    # Write replicate groups informations
    replicate_groups = dict()
    FH_replicate_log = open(replicate_log)
    for line in FH_replicate_log:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            replicate_groups[line[0]] = { 
            'Replicates' : line[1]
            }

    # Global after filters
    out_biom = BiomIO.from_json( output_biom )
    for observation_name in out_biom.get_observations_names():
        global_results['nb_clstr_kept'] += 1
        global_results['nb_seq_kept'] += out_biom.get_observation_count( observation_name )
    for sample_name in out_biom.get_samples_names():
        samples_results[sample_name]['kept'] = sum( 1 for x in out_biom.get_observations_by_sample(sample_name) )
    del out_biom

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "cluster_filters_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###PORCESSED_FILTERS###" in line:
            line = line.replace( "###PORCESSED_FILTERS###", json.dumps([filter for filter in discards]) )
        elif "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
        elif "###REPLICATE_GROUPS###" in line:
            line = line.replace( "###REPLICATE_GROUPS###", json.dumps(replicate_groups) )
        elif "###FILTERS_RESULTS###" in line:
            line = line.replace( "###FILTERS_RESULTS###", json.dumps(list(filters_results.values())) )
        elif "###DATA_SAMPLE###" in line:
            line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
        elif "###CLUSTERS_SIZES###" in line:
            line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
        elif "###DATA_COUNTS###" in line:
            line = line.replace( "###DATA_COUNTS###", json.dumps(counts) )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()

# def ratioParameter( arg_value ):
#     """
#     @summary: Argparse type for ratio (float between 0 and 1).
#     """
#     float_arg_value = None
#     try:
#         float_arg_value = float(arg_value)
#         if float_arg_value < 0.0 or float_arg_value > 1.0:
#             raise_exception( argparse.ArgumentTypeError("must be between 0.0 and 1.0."))
#     except:
#         raise_exception( argparse.ArgumentTypeError("must be between 0.0 and 1.0."))
#     return float_arg_value

def minAbundParameter( arg_value ):
    """
    @summary: Argparse type for min-abundance parameter.
    """
    str_value = arg_value.replace(",",".") 
    cleaned_value = float(str_value) if "." in str_value or "e-" in str_value else int(str_value)
    return cleaned_value



def process( args ):
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )

    try:
        discards = dict() # by filter the discard file path

        if args.min_sample_presence is not None:
            label = "Present in less than " + str(args.min_sample_presence) + " samples"
            discards[label] = tmpFiles.add( "min_sample_presence" )
            excluded_obs_on_samplePresence( args.input_biom, args.min_sample_presence, discards[label] )

        replicate_groups_log = tmpFiles.add( "replicate_groups.txt" )
        if args.min_replicate_presence is None:
            FH_log = Logger(replicate_groups_log)
            FH_log.write('No replicate groups defined\n')
        elif args.min_replicate_presence is not None and args.replicate_file is not None:
            label = "Present in less than " + str(args.min_replicate_presence*100) + "%  of replicates of all replicate groups."
            discards[label] = tmpFiles.add( "min_replicate_presence")
            excluded_obs_on_replicatePresence( args.input_biom, args.replicate_file, args.min_replicate_presence, replicate_groups_log, discards[label])

        if args.min_abundance is not None:
            
            if type(args.min_abundance) == float:
                biom = BiomIO.from_json( args.input_biom )
                min_nb_seq = int(biom.get_total_count() * args.min_abundance) + 1
                label = "Abundance < " + str(args.min_abundance*100) + "% (i.e " + str(min_nb_seq) + " sequences )"
            else:
                label = "Abundance < " + str(args.min_abundance)
            discards[label] = tmpFiles.add( "min_abundance" )
            excluded_obs_on_abundance( args.input_biom, args.min_abundance, discards[label] )

        if args.contaminant is not None:
            label = "Present in databank of contaminants"
            discards[label] = tmpFiles.add( "contaminant" )
            cleaned_seq = tmpFiles.add( "cleaned_sequences.fasta" )
            cleaned_biom = tmpFiles.add( "cleaned_abundance.biom" )
            contaminated_seq = tmpFiles.add( "contaminated_sequences.fasta" )
            remove_log = tmpFiles.add( "clean.log" )
            cmd = RemoveConta(args.contaminant, args.input_fasta, args.input_biom, cleaned_seq, cleaned_biom, contaminated_seq, remove_log, args.nb_cpus, args.debug)
            cmd.submit( args.log_file )
            cmd.write_conta_list( discards[label] )

        # Removes observations
        clusters_discarded = uniq_from_files_lists( [discards[filter] for filter in discards] )
        remove_observations( clusters_discarded, args.input_biom, args.output_biom )

        # Selects the N most abundant
        if args.nb_biggest_clusters:
            label = "Not in the " + str(args.nb_biggest_clusters) + " biggest"
            discards[label] = tmpFiles.add( "nb_biggest_clusters" )
            excluded_obs_on_nBiggest( args.output_biom, args.nb_biggest_clusters, discards[label] )
            # Removes observations
            remove_observations( uniq_from_files_lists([discards[label]]), args.output_biom, args.output_biom )

        # Writes outputs
        update_fasta_log = tmpFiles.add( "update_fasta_log.txt" )
        UpdateFasta( args.output_biom, args.input_fasta, args.output_fasta, update_fasta_log ).submit( args.log_file )
        write_exclusion( discards, args.excluded )
        
        depth_file = tmpFiles.add( "depths.tsv" )
        Depths(args.output_biom, depth_file).submit( args.log_file )
        
        write_summary( args.summary, args.input_biom, args.output_biom, replicate_groups_log, discards, depth_file )

    finally:
        if not args.debug : 
            tmpFiles.deleteAll()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description='Filters an abundance file')
    parser.add_argument('-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]")
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    #     Filters
    group_filter = parser.add_argument_group( 'Filters' )
    group_filter.add_argument( '--nb-biggest-clusters', type=int, default=None, required=False, help="Number of most abundant clusters you want to keep.") 
    group_filter.add_argument( '-s', '--min-sample-presence', type=int, help="Keep cluster present in at least this number of samples.") 
    group_filter.add_argument( '-r', '--min-replicate-presence', type=minAbundParameter, default=None, help="Keep cluster present in at least this proportion of replicates in at least one group (please indicate a proportion between 0 and 1). Replicates must be defined with --replicate_file REPLICATE FILE")
    group_filter.add_argument( '--replicate_file', help='Replicate file must be specified if --min-replicate-presence is set. First column of the file must indicate the sample name, and the second column the group name of this replicate. Exemple: TEM1_L0001_R   Temoin.')
    group_filter.add_argument( '-a', '--min-abundance', type=minAbundParameter, default=None, required=False, help="Minimum percentage/number of sequences, comparing to the total number of sequences, of a cluster (between 0 and 1 if percentage desired)." )
    # group_filter.add_argument( '--abundance-by-sample', type=bool, default=False, action='store_true', help="Abundance threshold is applied by default on the total abundance of cluster. Activate this option if you want to applied the threshold on sample abundances (if float, each cluster must be present in a " )
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help="The input BIOM file. (format: BIOM)")
    group_input.add_argument('--input-fasta', required=True, help="The input FASTA file. (format: FASTA)")
    group_input.add_argument('--contaminant', default=None, help="Use this databank to filter sequence before affiliation. (format: FASTA)")
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default="cluster_filters_abundance.biom", help="The BIOM file output. (format: BIOM) [Default: %(default)s]")
    group_output.add_argument('--output-fasta', default="cluster_filters.fasta", help="The FASTA output file. (format: FASTA) [Default: %(default)s]")
    group_output.add_argument('--summary', default="cluster_filters.html", help="The HTML file containing the graphs. [Default: %(default)s]")
    group_output.add_argument('--excluded', default="cluster_filters_excluded.tsv", help="The TSV file that summarizes all the clusters discarded. (format: TSV) [Default: %(default)s]")
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")

    if args.nb_biggest_clusters is None and args.min_sample_presence is None and args.min_replicate_presence is None and args.min_abundance is None and args.contaminant is None:
        raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : At least one filter must be set to run " + os.path.basename(sys.argv[0]) + "\n\n"))
    if not args.min_abundance is None and (args.min_abundance <= 0 or (type(args.min_abundance) == float and args.min_abundance >= 1.0 ) ):
        raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : If filtering on abundance, you must indicate a positive threshold and if percentage abundance threshold must be smaller than 1.0. \n\n" ))
    if not args.min_replicate_presence is None and args.replicate_file is None or not args.replicate_file is None and args.min_replicate_presence is None:
        raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : --min-replicate-presence and --replicate_file must be both specified.\n\n" ))
    elif not args.min_replicate_presence is None and (args.min_replicate_presence < 0 or args.min_replicate_presence > 1):
        raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : --min_replicate_presence must be a proportion (between 0 and 1).\n\n "))
    # Process
    process( args )
