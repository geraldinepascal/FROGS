#!/usr/bin/env python3

__author__ = 'Frédéric Escudié - Genotoul/MIAT'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
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


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def process( in_biom, out_biom, out_metadata ):
    ordered_blast_keys = ["taxonomy", "subject", "evalue", "perc_identity", "perc_query_coverage", "aln_length"] # Keys in blast_affiliations metadata
    taxonomy_depth = 0
    unclassified_observations = list()

    FH_metadata = open( out_metadata, "wt" )
    FH_metadata.write( "ASVID\t" + "\t".join([item for item in ordered_blast_keys]) + "\n" )
    biom = BiomIO.from_json( in_biom )
    for observation in biom.get_observations():
        for metadata_key in list(observation["metadata"].keys()):
            # Extract blast_affiliations metadata in metadata_file
            if metadata_key == "blast_affiliations": 
                if observation["metadata"][metadata_key] is not None:
                    for current_affi in observation["metadata"][metadata_key]:
                        if isinstance(current_affi["taxonomy"], list) or isinstance(current_affi["taxonomy"], tuple):
                            current_affi["taxonomy"] = ";".join( current_affi["taxonomy"] )
                        FH_metadata.write( observation["id"] + "\t" + "\t".join([str(current_affi[item]) for item in ordered_blast_keys]) + "\n" )
                del observation["metadata"][metadata_key]
            # All list are transformed in string (including blast_taxonomy)
            elif observation["metadata"][metadata_key] is not None: 
                if isinstance(observation["metadata"][metadata_key], list) or isinstance(observation["metadata"][metadata_key], tuple):
                    observation["metadata"][metadata_key] = ";".join( map(str, observation["metadata"][metadata_key]) )
        if "blast_taxonomy" in observation["metadata"]:
            if observation["metadata"]["blast_taxonomy"] is None or len(observation["metadata"]["blast_taxonomy"])==0:
                unclassified_observations.append( observation["id"] )
                observation["metadata"]["taxonomy"] = list()
            else:
                taxonomy_depth = len(observation["metadata"]["blast_taxonomy"].split(";"))
                observation["metadata"]["taxonomy"] = observation["metadata"]["blast_taxonomy"].split(";")
    # Add "Unclassified" ranks in unclassified observations
    if taxonomy_depth > 0:
        for observation_id in unclassified_observations:
            observation_metadata = biom.get_observation_metadata(observation_id)
            observation_metadata["taxonomy"] = ["Unclassified"] * taxonomy_depth
    BiomIO.write( out_biom, biom )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='The detailed FROGS blast affiliations can trigger problem with tools like Qiime. This script extracts the problematic metadata in a second file and writes a BIOM usable in every tool using BIOM.' )
    parser.add_argument('--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help="The abundance file (format: BIOM)." )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default='abundance.std.biom', help='The fully compatible abundance file (format: BIOM). [Default: %(default)s]' )
    group_output.add_argument('--output-metadata', default='blast_informations.std.tsv', help='The blast affiliations informations (format: TSV). [Default: %(default)s]' )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several information on executed commands. [Default: stdout]' )
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Process
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    process( args.input_biom, args.output_biom, args.output_metadata )
