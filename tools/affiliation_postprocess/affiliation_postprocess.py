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

__author__ = 'Maria Bernard INRA - SIGENAE'
__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
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
from frogsSequenceIO import *

###################################################################################################################
###                                           CLASS                                                             ###
###################################################################################################################

class SelectInclusiv(Cmd):
    """
    @summary : Select smallest amplicon reference among multiaffiliations
    """
    def __init__(self, reference, biom_in, biom_out, log ):
        """
        @param reference : [str] Path targeted amplicon reference fasta file
        @param biom_in : [str] Path to biom file to work on
        @param biom_out : [str] Path to biom file to output
        @param log : [str] Path to log file
        """
        Cmd.__init__( self,
                      'select_inclusive_amplicon.py',
                      'Select smallest reference as affiliation.',
                      ' --ITS-reference ' + reference + ' -b ' + biom_in + ' -o ' + biom_out + " --log-file " + log,
                      '--version' )
        self.inclusiv_log = log

    def parser(self, log_file):
        FH_in = open(self.inclusiv_log)
        FH_log = Logger( log_file )
        for line in FH_in:
            FH_log.write( '\t'+line )
        FH_log.write('\n')
        FH_log.close()
        FH_in.close()


    def get_version(self):   
        return Cmd.get_version(self, 'stderr').strip()

class AggregateOTU(Cmd) : 
    """
    @summary : aggregate OTU that share at least one taxonomic affiliation
    """
    def __init__(self, input_biom, input_fasta, identity, coverage, ignored_list, output_biom, output_compo, output_fasta, log):
        """
        @param input_biom : [str] Path to biom file to treat
        @param input_fasta : [str] Path to fasta file to treat
        @param identity : [float] min percentage identity of taxonomic affiliations to merge OTU
        @param coverage : [float] min percentage coverage of taxonomic affiliations to merge OTU
        @param ignored_list : [str] list of taxon ignored for OTU aggregation
        @param output_biom : [str] Path to the resulting biom file
        @param output_compo : [str] Path to the resulting TSV file, containing aggregated OTU composition
        @param output_fasta : [str] Path to the resulting fasta file
        @param log : [str] Path to excution log file
        """
        opt = ""
        if ignored_list:
            opt += ' --taxon-ignored ' + ' '.join(['\"' + t + '\"' for t in ignored_list]) + ''
        Cmd.__init__( self,
                      'aggregate_affiliated_otus.py',
                      'Aggregate OTU that share taxonomic affiliation with at least I% identity and C% coverage',
                      ' -b ' + input_biom + ' -f ' + input_fasta + opt + ' --identity ' + str(identity) + ' --coverage ' + str(coverage) \
                      + ' --output-biom ' + output_biom +' --output-compo ' +  output_compo + ' --output-fasta ' + output_fasta + " --log-file " + log,
                      '--version' )

        self.aggregate_log = log

    def parser(self, log_file):
        FH_in = open(self.aggregate_log)
        FH_log = Logger( log_file )
        for line in FH_in:
            FH_log.write( '\t'+line )
        FH_log.write('\n')
        FH_log.close()
        FH_in.close()

    def get_version(self):   
        return Cmd.get_version(self, 'stderr').strip()

###################################################################################################################
###                                           FUNCTIONS                                                         ###
###################################################################################################################

def process(params):
    """
    @summary : refine affiliation by:
        - keep the smallest ITS amplicon refence for OTU multiaffiliated, these OTU are most likely inclusive ITS
        - agggregate OTU that share the same affiliation based on threshold on %coverage and %id
    """
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )
    if params.reference:
        smallest_its_biom = tmpFiles.add(os.path.basename(args.input_biom + "_resolve_inclusiv_its.biom"))
        smallest_its_log = tmpFiles.add(os.path.basename(args.input_biom + "_resolve_inclusiv_its.log"))
    aggregate_tmp_log = tmpFiles.add(os.path.basename(args.input_biom + "_aggregation.log"))
    input_biom = params.input_biom
    try:
        Logger.static_write(params.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
        # inclusive ITS case
        if params.reference:
            select_inclusiv_cmd = SelectInclusiv(params.reference, input_biom, smallest_its_biom, smallest_its_log)
            select_inclusiv_cmd.submit(params.log_file)
            input_biom = smallest_its_biom

        # aggregate OTU
        aggregate_cmd = AggregateOTU(input_biom, params.input_fasta, params.identity, params.coverage, params.taxon_ignored, params.output_biom, params.output_compo, params.output_fasta, aggregate_tmp_log)
        aggregate_cmd.submit(params.log_file)
    finally:
        if not args.debug:
            tmpFiles.deleteAll()

###################################################################################################################
###                                           MAIN                                                              ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Refine affiliations, to manage amplicon included in other sequence, and to deal with surnumerary OTU (OTU with same affiliations.")
    parser.add_argument( '-i', '--identity', default=99.0, help="Min percentage identity to agggregate OTU. [Default: %(default)s]")
    parser.add_argument( '-c', '--coverage', default=99.0, help="Min percentage coverage to agggregate OTU. [Default: %(default)s]")
    parser.add_argument( '-t', '--taxon-ignored', type=str, nargs='*', help="Taxon list to ignore when OTUs agggregation")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)

    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-b', '--input-biom', required=True, help='Abundance table with affiliations metadata from the affiliation_OTU program (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help='OTU seed sequence file (format: FASTA).')
    group_input.add_argument('-r', '--reference', required=False, help='amplicon reference file, to resolve inclusive amplicon affiliations (format: FASTA)')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', default='affiliation_postprocess_abundance.biom', help='BIOM file whith refind affiliation annotations. (format: BIOM) [Default: %(default)s]')
    group_output.add_argument('--output-compo', default='affiliation_postprocess_otu_composition.tsv', help='Aggregated OTU composition (format: TSV) [Default: %(default)s]')
    group_output.add_argument('--output-fasta', default='affiliation_postprocess_OTU.fasta', help='Updated OTU FASTA file (format: FASTA) [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)
