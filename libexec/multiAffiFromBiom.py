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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import BiomIO
from frogsUtils import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def multiAffiFromBiom( biom, output_tsv ):
    """
    @summary: Extracts multi-affiliations data from a FROGS BIOM file.
    @param input_biom: [Biom] The BIOM object.
    @param output_tsv: [str] Path to the output file (format : TSV).
    """
    out_fh = open( output_tsv, "wt" )
#     out_fh.write( "#" + "\t".join(["OTU", "Subject_taxonomy", "Blast_subject", "Prct_identity", "Prct_query_coverage", "e-value", "Alignment_length"]) + "\n" )
    out_fh.write( "#" + "\t".join(["observation_name", "blast_taxonomy", " blast_subject", "blast_perc_identity", "blast_perc_query_coverage", "blast_evalue", "blast_aln_length"]) + "\n" )
    for current_observation in biom.get_observations():
        if issubclass(current_observation['metadata']["blast_affiliations"].__class__, list) and len(current_observation["metadata"]["blast_affiliations"]) > 1:
            for current_aln in current_observation["metadata"]["blast_affiliations"]:
                taxonomy = current_aln["taxonomy"]
                if isinstance(taxonomy, list) or isinstance(taxonomy, tuple):
                    taxonomy = ";".join( current_aln["taxonomy"] )
                line_parts = [
                    current_observation["id"],
                    taxonomy,
                    str(current_aln["subject"]),
                    str(current_aln["perc_identity"]),
                    str(current_aln["perc_query_coverage"]),
                    str(current_aln["evalue"]),
                    str(current_aln["aln_length"])
                ]
                out_fh.write( "\t".join(line_parts) + "\n" )
    out_fh.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Extracts multi-affiliations data from a FROGS BIOM file.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file', required=True, help='Path to the abundance file (format: BIOM).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='Path to the output file (format: TSV).')
    args = parser.parse_args()

    # Process
    biom = BiomIO.from_json( args.input_file )
    if not biom.has_observation_metadata( 'blast_affiliations' ):
        raise_exception( Exception( "\n\n#ERROR : '" + args.input_file  + "' cannot be used in " + sys.argv[0] + ": this file does not contain 'blast_affiliations'.\n\n" ))
    multiAffiFromBiom( biom, args.output_file )
