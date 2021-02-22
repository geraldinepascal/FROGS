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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse and Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2.1'
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
# COMMAND LINES
#
##################################################################################################################################################
class Biom2tsv(Cmd):
    """
    @summary: Converts BIOM file to TSV file.
    @note: taxonomyRDP seedID seedSequence blastSubject blastEvalue blastLength blastPercentCoverage blastPercentIdentity blastTaxonomy OTUname SommeCount sample_count
    """
    def __init__( self, out_tsv, in_biom, headerOnly, in_fasta=None ):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv: [str] Path to output TSV file.
        """
        # Sequence file option
        sequence_file_opt = "" if in_fasta is None else " --input-fasta " + in_fasta

        # Check the metadata
        biom = BiomIO.from_json( in_biom )
        obs = biom.rows[0]
        conversion_tags = ""
        if biom.has_observation_metadata( 'comment' ) :
            conversion_tags += "'comment' "
        if biom.has_observation_metadata( 'status' ) :
            conversion_tags += "'status' "
        if biom.has_observation_metadata( 'rdp_taxonomy' ) and biom.has_observation_metadata( 'rdp_bootstrap' ):
            conversion_tags += "'@rdp_tax_and_bootstrap' "
        if biom.has_observation_metadata( 'blast_taxonomy' ):
            conversion_tags += "'blast_taxonomy' "
        if biom.has_observation_metadata( 'blast_affiliations' ):
            conversion_tags += "'@blast_subject' "
            conversion_tags += "'@blast_perc_identity' "
            conversion_tags += "'@blast_perc_query_coverage' "
            conversion_tags += "'@blast_evalue' "
            conversion_tags += "'@blast_aln_length' "
        if biom.has_observation_metadata( 'seed_id' ):
            conversion_tags += "'seed_id' "
        if in_fasta is not None:
            conversion_tags += "'@seed_sequence' "

        frogs_metadata = ["comment", 'status' , "rdp_taxonomy", "rdp_bootstrap","blast_taxonomy","blast_affiliations","seed_id"]
        if biom.get_observation_metadata(obs["id"]) != None:
            for metadata in biom.get_observation_metadata(obs["id"]):
                if metadata not in frogs_metadata : 
                    conversion_tags += "'"+metadata+"' "
                
        conversion_tags += "'@observation_name' '@observation_sum' '@sample_count'"

        if headerOnly:
            header_list = conversion_tags.replace("'","").replace("@","").split()[0:-1]
            for sample_name in biom.get_samples_names():
                header_list.append(sample_name)

            Cmd.__init__( self,
                          'echo',
                          'Print biom metadata as TSV header file',
                          '\"#' + "\\t".join(header_list) + " \" > " + out_tsv )
        else:
            # Set command
            Cmd.__init__( self,
                          'biom2tsv.py',
                          'Converts a BIOM file in TSV file.',
                          "--input-file " + in_biom + sequence_file_opt + " --output-file " + out_tsv + " --fields " + conversion_tags,
                          '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()                          


class Biom2multiAffi(Cmd):
    """
    @summary: Extracts multi-affiliations from a FROGS BIOM file.
    """
    def __init__( self, out_tsv, in_biom, headerOnly ):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv: [str] Path to output TSV file.
        """
        if headerOnly:
            header_list = ["observation_name", "blast_taxonomy", " blast_subject", "blast_perc_identity", "blast_perc_query_coverage", "blast_evalue", "blast_aln_length"]
            Cmd.__init__( self,
                          'echo',
                          'Print biom blast multiAffiliation as TSV header file',
                          '\"#' + "\\t".join(header_list) + " \" > " + out_tsv )
        else:
            Cmd.__init__( self,
                          'multiAffiFromBiom.py',
                          'Extracts multi-affiliations data from a FROGS BIOM file.',
                          '--input-file ' + in_biom + ' --output-file ' + out_tsv,
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
    parser = argparse.ArgumentParser( description='Converts a BIOM file in TSV file.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '--header', default=False, action='store_true', help="Print header only" )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-b', '--input-biom', required=True, help="The abundance file (format: BIOM)." )
    group_input.add_argument( '-f', '--input-fasta', help='The sequences file (format: FASTA). If you use this option the sequences will be add in TSV.' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-t', '--output-tsv', default='abundance.tsv', help='This output file will contain the abundance and metadata (format: TSV). [Default: %(default)s]' )
    group_output.add_argument( '-m', '--output-multi-affi', default='multihits.tsv', help='This output file will contain information about multiple alignements (format: TSV). Use this option only if your affiliation has been produced by FROGS. [Default: %(default)s]' )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands.' )
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Process
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    Biom2tsv( args.output_tsv, args.input_biom, args.header, args.input_fasta ).submit( args.log_file )
    
    biom = BiomIO.from_json( args.input_biom )
    if biom.has_metadata("blast_affiliations"):
        Biom2multiAffi( args.output_multi_affi, args.input_biom, args.header ).submit( args.log_file )
