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

__author__ = 'Olivier Rue - Migale MaIAGE Jouy-en-Josas - Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
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
from frogsBiom import BiomIO
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class ITSx(Cmd):
    """
    @summary: Use ITSx to identifies ITS sequences and extracts the ITS region
    """
    def __init__(self, in_fasta, in_biom, target, out_fasta, out_count, out_removed, log_file, param):
        """
        @param in_fasta: [str] Path to the fasta to process.
        @param in_count: [str] Path to the associated count file to update.
        @param target : [str] Either ITS1 or ITS2
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
        Cmd.__init__(self,
            'parallelITSx.py',
            'identifies ITS sequences and extracts the ITS region',
            ' -f ' + in_fasta + ' -b ' + in_biom + options + ' --its '+ target + ' -o ' + out_fasta + ' -m ' + out_removed + ' -a ' + out_count + ' --log-file ' + log_file,
            '--version'
            )
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

def write_summary( summary_file, input_biom, output_biom ):
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

    # print samples_results

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "itsx_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
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
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '--region', type=str, required=True,  choices=['ITS1','ITS2'], help='Which fungal ITS region is targeted: either ITS1 or ITS2' )
    group_input.add_argument( '--check-its-only', action='store_true', default=False, help='Check only if sequences seem to be an ITS' )
    group_input.add_argument( '-f', '--input-fasta', required=True, help='The cluster sequences (format: FASTA).' )
    group_exclusion_abundance = group_input.add_mutually_exclusive_group()
    group_exclusion_abundance.add_argument( '-b', '--input-biom', help='The abundance file for clusters by sample (format: BIOM).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-n', '--out-fasta', default='itsx.fasta', help='sequences file out from ITSx (format: FASTA). [Default: %(default)s]')
    group_output.add_argument( '-a', '--out-abundance', default="itsx_abundance.biom", help='Abundance file without chimera (format: BIOM ). [Default: %(default)s]')
    group_output.add_argument( '-m', '--out-removed', default='itsx_removed.fasta', help='sequences file removed (format: FASTA). [Default: %(default)s]')
    group_output.add_argument( '--summary', default="itsx.html", help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.out_fasta)[0] )
    
    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        log_itsx = tmpFiles.add("ITSx.log")
        
        ITSx(args.input_fasta, args.input_biom, args.region, args.out_fasta, args.out_abundance, args.out_removed, log_itsx, args ).submit( args.log_file )
        write_summary( args.summary, args.input_biom, args.out_abundance)
        
        # Append independant log files
        log_append_files( args.log_file, [log_itsx] )
        
    # Remove temporary files
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
