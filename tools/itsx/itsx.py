#!/usr/bin/env python3

__author__ = 'Olivier Rué - Migale/MaIAGE & Maria Bernard - SIGENAE/GABI'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import json
import argparse

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

class ITSx(Cmd):
    """
    @summary: Use ITSx to identifies ITS sequences and extracts the ITS region
    """
    def __init__(self, in_fasta, in_biom, organism_groups, out_fasta, out_count, out_removed, log_file, param):
        """
        @param in_fasta: [str] Path to the fasta to process.
        @param in_count: [str] Path to the associated count file to update.
        @param target : [str] Either ITS1 or ITS2
        @param organism_groups: [list] organism groups to scan (default only Fungi)
        @param out_fasta: [str] Path to the processed fasta.
        @param out_count: [str] Path to the updated count file.
        @param log_file: [str] Path to the updated count file.
        @param param: [Namespace] The 'param.nb_cpus'.
        """
        options=''
        if param.nb_cpus > 1:
            options = " --nb-cpus " + str(param.nb_cpus)
        if param.debug:
            options += " --debug "
        if param.check_its_only:
            options += " --check-its-only"
        else :
            options += " --its "+param.region
        Cmd.__init__(self,
            'parallelITSx.py',
            'identifies ITS sequences and extracts the ITS region',
            ' -f ' + in_fasta + ' -b ' + in_biom + options + ' --organism-groups ' + ' '.join(organism_groups) + ' -o ' + out_fasta + ' -m ' + out_removed + ' -a ' + out_count + ' --log-file ' + log_file,
            '--version')
        self.program_log = log_file


    def parser(self, log_file):
        """
        @summary: Parse the command results to add information in log_file.
        @log_file: [str] Path to the ITSx log file.
        """
        # Parse output
        FH_log_ITSX = open( self.program_log )
        count_ITSx = dict()
        kept = ""
        for line in FH_log_ITSX:
            if line.startswith('nb') :
                [key,value]=line.strip().split(':')
                if key in count_ITSx:
                    count_ITSx[key]+= int(value)
                else:
                    count_ITSx[key]= int(value)
                if "kept" in key:
                    kept = key
        FH_log_ITSX.close()
        # Write result
        FH_log = Logger( log_file )
        FH_log.write( 'Results:\n' )
        #FH_log.write( '\t' + kept + ' : ' + str(count_ITSx[kept]) + '\n')
        #for key in count_ITSx :
        #    if key != kept:
        #        FH_log.write( '\t' + key + ' : ' + str(count_ITSx[key]) + '\n')

        FH_log.close()

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').strip()


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
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

def write_summary( summary_file, input_biom, output_biom, depth_file ):
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
    
    
    
    ## Retrieve IDs to remove
    # Get initial observation names
    in_biom = BiomIO.from_json(input_biom)
    all_clusters_ids = list()
    for observation_name in in_biom.get_observations_counts():
        all_clusters_ids.append(str(observation_name[0]))
    # Get kept observation names
    out_biom = BiomIO.from_json( output_biom )
    kept_clusters_ids = list()
    for observation_name in out_biom.get_observations_counts():
        kept_clusters_ids.append(str(observation_name[0]))
    # Difference between inital and kept
    lost_clusters = [x for x in all_clusters_ids if x not in kept_clusters_ids]
    
    
    
    # Global before filters
    in_biom = BiomIO.from_json( input_biom )
    for observation_name in in_biom.get_observations_names():
        global_results['nb_clstr_ini'] += 1
        global_results['nb_seq_ini'] += in_biom.get_observation_count( observation_name )
    for sample_name in in_biom.get_samples_names():
        samples_results[sample_name] = {
            'initial': sum( 1 for x in in_biom.get_observations_by_sample(sample_name) ),
            'filtered': dict(),
            'kept': 0,
            'initial_ab': sum ( in_biom.get_count( observation['id'], sample_name ) for observation in in_biom.get_observations_by_sample(sample_name)),
            'kept_ab':0
        }
    
    # Global after filters
    out_biom = BiomIO.from_json( output_biom )
    output_biom_test = BiomIO.from_json( output_biom )
    for observation_name in out_biom.get_observations_names():
        global_results['nb_clstr_kept'] += 1
        global_results['nb_seq_kept'] += out_biom.get_observation_count( observation_name )
    for sample_name in out_biom.get_samples_names():
        samples_results[sample_name]['kept'] = sum( 1 for x in out_biom.get_observations_by_sample(sample_name) )
        samples_results[sample_name]['kept_ab'] = sum ( output_biom_test.get_count( obs['id'], sample_name ) for obs in output_biom_test.get_observations_by_sample(sample_name))
    del out_biom

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

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "itsx_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###DATA_SAMPLE###" in line:
            line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
        elif "###CLUSTERS_SIZES###" in line:
            line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
        elif "###DATA_COUNTS###" in line:
            line = line.replace( "###DATA_COUNTS###", json.dumps(counts) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
        elif "###FROGS_VERSION###" in line:
            line = line.replace( "###FROGS_VERSION###", "\""+str(__version__)+"\"" )
        elif "###FROGS_TOOL###" in line:
            line = line.replace( "###FROGS_TOOL###", "\""+ os.path.basename(__file__)+"\"" )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()




##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Removes PCR chimera.'
    )
    parser.add_argument('--version', action='version', version=__version__ )
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )
    parser.add_argument('--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    
    # Params
    group_params = parser.add_argument_group( 'Parameters' )
    group_params.add_argument( '--organism-groups', type=str, nargs="+", default=['F'], help='Reduce ITSx scan to specified organim groups. [Default: %(default)s , which means Fungi only]')
    group_exclusion_params = group_params.add_mutually_exclusive_group(required=True)
    group_exclusion_params.add_argument( '--region', type=str, default=None, choices=['ITS1','ITS2'], help='Which fungal ITS region is targeted and trimmed: either ITS1 or ITS2. (mutually exclusive with --check-its-only) [Default: %(default)s]' )
    group_exclusion_params.add_argument( '--check-its-only', action='store_true', default=False, help='Check only if sequences seem to be an ITS (mutually exclusive with --region) [Default: %(default)s]' )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-fasta', required=True, help='The cluster sequences (format: FASTA).' )
    group_input.add_argument('--input-biom', help='The abundance file for clusters by sample (format: BIOM).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-fasta', default='itsx.fasta', help='sequences file out from ITSx (format: FASTA). [Default: %(default)s]')
    group_output.add_argument('--output-biom', default="itsx_abundance.biom", help='Abundance file without chimera (format: BIOM ). [Default: %(default)s]')
    group_output.add_argument('--output-removed-sequences', default='itsx_removed.fasta', help='sequences file removed (format: FASTA). [Default: %(default)s]')
    group_output.add_argument('--html', default="itsx.html", help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.output_fasta)[0] )
    
    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        log_itsx = tmpFiles.add("ITSx.log")
        
        ITSx(args.input_fasta, args.input_biom, args.organism_groups, args.output_fasta, args.output_biom, args.output_removed_sequences, log_itsx, args ).submit( args.log_file )
        
        depth_file = tmpFiles.add( "depths.tsv" )
        Depths(args.output_biom, depth_file).submit( args.log_file )
        
        write_summary( args.html, args.input_biom, args.output_biom, depth_file)
        
        # Append independant log files
        log_append_files( args.log_file, [log_itsx] )
        
    # Remove temporary files
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
