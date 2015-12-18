#!/usr/bin/env python2.7
#
# Copyright (C) 2015 INRA
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

__author__ = 'Maria Bernard - Sigenae AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.3.3'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import sys
import argparse
from biom import *
from sequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def uniq( metadata_list, evaluated_tag, ambiguity_value ):
    value = None
    for metadata_elt in metadata_list:
        if value is None:
            value = metadata_elt[evaluated_tag]
        elif value != metadata_elt[evaluated_tag]:
            value = ambiguity_value
    return value

def observation_line_parts( observation, count_by_sample, fields, list_separator ):
    no_affiliation_str = "no data"
    line = list()
    for current_field in fields:
        if current_field == '@observation_name':
            line.append( str(observation['id']) )
        elif current_field == '@sample_count':
            line.append( "\t".join(map(str, count_by_sample)) )
        elif current_field == '@observation_sum':
            line.append( str(sum(count_by_sample)) )
        elif current_field == "@rdp_tax_and_bootstrap":
                rdp_and_bootstrap = ""
                rdp_taxonomy = observation['metadata']["rdp_taxonomy"]
                rdp_bootstrap = observation['metadata']["rdp_bootstrap"]
                for i, tax in enumerate(rdp_taxonomy):
                    rdp_and_bootstrap += tax + ";(" + str(rdp_bootstrap[i]) + ");" # tax1;(boots1);tax2;(boots2);
                line.append(str(rdp_and_bootstrap))
        elif current_field == "@blast_perc_identity":
            if len(observation['metadata']["blast_affiliations"]) > 0:
                line.append( str(uniq(observation['metadata']["blast_affiliations"], "perc_identity", "multi-identity")) )
            else:
                line.append( no_affiliation_str )
        elif current_field == "@blast_perc_query_coverage":
            if len(observation['metadata']["blast_affiliations"]) > 0:
                line.append( str(uniq(observation['metadata']["blast_affiliations"], "perc_query_coverage", "multi-coverage")) )
            else:
                line.append( no_affiliation_str )
        elif current_field == "@blast_evalue":
            if len(observation['metadata']["blast_affiliations"]) > 0:
                line.append( str(uniq(observation['metadata']["blast_affiliations"], "evalue", "multi-evalue")) )
            else:
                line.append( no_affiliation_str )
        elif current_field == "@blast_subject":
            if len(observation['metadata']["blast_affiliations"]) > 0:
                line.append( str(uniq(observation['metadata']["blast_affiliations"], "subject", "multi-subject")) )
            else:
                line.append( no_affiliation_str )
        elif current_field == "@blast_aln_length":
            if len(observation['metadata']["blast_affiliations"]) > 0:
                line.append( str(uniq(observation['metadata']["blast_affiliations"], "aln_length", "multi-alignment-lg")) )
            else:
                line.append( no_affiliation_str )
        else: #metadata
            if issubclass(observation['metadata'][current_field].__class__, list):
                line.append( list_separator.join(observation['metadata'][current_field]) )
            else:
                line.append( str(observation['metadata'][current_field]) )
    return line

def header_line_parts( fields, biom ):
    header_parts = list()
    for current_field in fields:
        if current_field == '@observation_name':
            header_parts.append( "observation_name" )
        elif current_field == '@sample_count':
            header_parts.append( "\t".join(biom.get_samples_names()) )
        elif current_field == '@observation_sum':
            header_parts.append( "observation_sum" )
        elif current_field == '@rdp_tax_and_bootstrap':
            header_parts.append( "rdp_tax_and_bootstrap" )
        elif current_field == "@blast_perc_identity":
            header_parts.append( "blast_perc_identity" )
        elif current_field == "@blast_perc_query_coverage":
            header_parts.append( "blast_perc_query_coverage" )
        elif current_field == "@blast_evalue":
            header_parts.append( "blast_evalue" )
        elif current_field == "@blast_subject":
            header_parts.append( "blast_subject" )
        elif current_field == "@blast_aln_length":
            header_parts.append( "blast_aln_length" )
        elif current_field == '@seed_sequence':
            header_parts.append( "seed_sequence" )
        else: #metadata
            header_parts.append( str(current_field) )
    return header_parts

def biom_to_tsv( input_biom, output_tsv, fields, list_separator ):
    """
    @summary: Convert BIOM file to TSV file.
    @param input_biom: [str] Path to the BIOM file.
    @param output_tsv: [str] Path to the output file (format : TSV).
    @param fields: [list] Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count' '@rdp_tax_and_bootstrap' . The others columns must be metadata title.
    @param list_separator: [str] Separator for complex metadata.
    """
    biom = BiomIO.from_json( input_biom )
    out_fh = open( output_tsv, "w" )
    # Header
    header_parts = header_line_parts( fields, biom )
    out_fh.write( "#" + "\t".join(header_parts) + "\n" )
    # Data
    for obs_idx, count_by_sample in enumerate(biom.to_count()):
        observation_parts = observation_line_parts( biom.rows[obs_idx], count_by_sample, fields, list_separator )
        out_fh.write( "\t".join(observation_parts) + "\n" )
    out_fh.close()

def biom_fasta_to_tsv( input_biom, input_fasta, output_tsv, fields, list_separator ):
    """
    @summary: Convert BIOM file to TSV file with sequence.
    @param input_biom: [str] Path to the BIOM file.
    @param input_fasta: [str] Path to the sequences of the observations.
    @param output_tsv: [str] Path to the output file (format : TSV).
    @param fields: [list] Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count', '@rdp_tax_and_bootstrap', '@seed_sequence'. The others columns must be metadata title.
    @param list_separator: [str] Separator for complex metadata.
    """
    biom = BiomIO.from_json( input_biom )
    out_fh = open( output_tsv, "w" )
    sequence_idx = fields.index("@seed_sequence")
    # Header
    header_parts = header_line_parts( fields, biom )
    out_fh.write( "#" + "\t".join(header_parts) + "\n" )
    # Data
    fields_without_seq = fields
    del fields_without_seq[sequence_idx]
    FH_in = FastaIO( input_fasta )
    for record in FH_in:
        obs_idx = biom.find_idx("observation", record.id)
        count_by_sample = biom.data.get_row_array(obs_idx)
        observation_parts = observation_line_parts( biom.rows[obs_idx], count_by_sample, fields_without_seq, list_separator )
        observation_parts.insert( sequence_idx, record.string )
        out_fh.write( "\t".join(observation_parts) + "\n" )
    out_fh.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Convert BIOM file to TSV file.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--fields', default=['@observation_name', '@observation_sum', '@sample_count'], nargs='+', help="Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count' '@rdp_tax_and_bootstrap', '@seed_sequence' . The others columns must be metadata title.")
    parser.add_argument( '-s', '--list-separator', default=';', help='Separator for complex metadata.')
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file', required=True, help='Path to the BIOM file.' )
    group_input.add_argument( '-a', '--input-fasta', default=None, required=False, help='Path to the FASTA file.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='Path to the output file (format : TSV).')
    args = parser.parse_args()

    # Process
    if args.input_fasta is not None:
        biom_fasta_to_tsv( args.input_file, args.input_fasta, args.output_file, args.fields, args.list_separator )
    else:
        biom_to_tsv( args.input_file, args.output_file, args.fields, args.list_separator )