#!/usr/bin/env python3

__author__ = 'Maria Bernard - SIGENAE/GABI'
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


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Tsv2biom(Cmd):
    """
    @summary: Converts TSV file to Biom file.
    @note: sample names are column name after observation_name followed by observation_sums. Everything before is considered as cluster metadata
    """
    def __init__( self, in_tsv, in_multi, out_biom, out_fasta=None ):
        """
        @param in_tsv: [str] Path to TSV file.
        @param in_multi: [str] Path to TSV multiaffiliation file.
        @param out_biom: [str] Path to output BIOM file.
        @param out_fasta: [str] Path to output FASTA file.
        """
        FH=open(in_tsv)
        header=FH.readline()
        FH.close()
        if header.startswith('#'):
            header=header[1:]
        header=header.strip().replace('"','').split()

        # Sequence file option
        if not out_fasta is None:
            sequence_file_opt = " --output-fasta " + out_fasta
        else:
            sequence_file_opt = ""

        if not in_multi is None:
            multihit_file_opt = " --input-multihits "+in_multi
        else:
            multihit_file_opt = ""

        # Check the metadata
        observation_name_index = header.index("observation_name")
        observation_sum_index = header.index("observation_sum")
        if (observation_sum_index - observation_name_index) != 1:
            raise_exception( Exception( "\n\nYou change the order of columns. TSV file must ended with observation_name, observation_sum, sample1, sample2 ... \n\n" ))
        samples_names = " ".join(header[observation_sum_index+1:])
        fields = " ".join(header[:observation_sum_index])
        # Set command
        Cmd.__init__( self,
                      'tsv2biom.py',
                      'Converts a TSV file in Biom.',
                      "--input-file "  + in_tsv + multihit_file_opt + " --output-file " + out_biom + sequence_file_opt + " --fields " + fields + " --samples-names " + samples_names,
                      '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()                      
        

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Converts a TSV file in BIOM file.' )
    parser.add_argument('--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-tsv', required=True, help='This input file contain the abundance and metadata (format: TSV).' )
    group_input.add_argument('--input-multi-affi', default=None, help='This input file will contain information about multiple alignements (format: TSV). Use this option only if your affiliation has been produced by FROGS. [Default: %(default)s]' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default='abundance.biom', help="The output abundance file (format: BIOM). [Default: %(default)s]" )
    group_output.add_argument('--output-fasta', default=None, help='The output sequences file (format: FASTA). If sequences exist in your input TSV with tag seed_sequence. [Default: %(default)s]' )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]' )
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Process
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    Tsv2biom( args.input_tsv, args.input_multi_affi, args.output_biom, args.output_fasta ).submit( args.log_file )
