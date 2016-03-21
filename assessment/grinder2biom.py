#!/usr/bin/env python2.7
#
# Copyright (C) 2016 INRA
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

__author__ = 'Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'beta'

import re
import argparse
from biom import Biom, BiomIO
from sequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def getCleanedTaxonomy( taxonomy ):
    """
    @summary: Returns the taxonomic ranks for the observation.
    @param observation_name: [str] The specified observation.
    @param taxonomy_key: [str] The metadata key for the taxonomy.
    @retrurn: [list] The list of taxonomic ranks.
    @note: Unfortunately some BIOM have a non-canonical format for store the taxonomy metadata. This method manages the below formats.
        - list or tuple:
          ["d:Bacteria", "Proteobacteria", "Epsilonproteobacteria", "Campylobacterales", "Helicobacteraceae", "Helicobacter"]
        - string:
           "Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter"
        - string ended by rank separator:
           "Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Helicobacter;"
        - string with bootstrap:
           "Bacteria(1.0000);Proteobacteria(0.9997);Epsilonproteobacteria(1.0000);Campylobacterales(1.0000);Helicobacteraceae(0.9898);Helicobacter(0.9912)"
        - string with bootstrap and ended by rank separator:
           "Bacteria(1.0000);Proteobacteria(0.9997);Epsilonproteobacteria(1.0000);Campylobacterales(1.0000);Helicobacteraceae(0.9898);Helicobacter(0.9912);"
    """
    cleaned_taxonomy = list()
    # Get taxonomy as a r/w list
    if isinstance(taxonomy, list) or isinstance(taxonomy, tuple): # Copy the list
        cleaned_taxonomy = [taxon.strip() for taxon in taxonomy]
    else: # Convert taxonomy in list
        cleaned_taxonomy = taxonomy
        if cleaned_taxonomy.strip().endswith(";"):
            cleaned_taxonomy = cleaned_taxonomy.strip()[:-1]
        if len(cleaned_taxonomy.split(";")) <= 3: # The tax separator is ","
            cleaned_taxonomy = cleaned_taxonomy.replace(",", ";")
        cleaned_taxonomy = [taxon.strip() for taxon in cleaned_taxonomy.split(";")]
    # Remove bootstrap information if its exist
    boostrap_regexp = re.compile("^(.+)\(\d+(\.\d+)?\)$")
    if len(cleaned_taxonomy) != 0 and boostrap_regexp.match(cleaned_taxonomy[0]) is not None: # Taxonomy contains bootstrap values
        for rank, taxon in enumerate(cleaned_taxonomy):
            matches = boostrap_regexp.search(taxon)
            cleaned_taxonomy[rank] = matches.group(1).strip()
    # Remove IDs
    cleaned_taxonomy = [taxon.split('[id:')[0].strip() for taxon in cleaned_taxonomy] # remove "[id: .....]"
    # Remove quotes
    for rank, taxon in enumerate(cleaned_taxonomy):
        cleaned_taxonomy[rank] = cleaned_taxonomy[rank].replace('\"', "").replace('"', "")
    # Remove root
    if cleaned_taxonomy[0].lower() == "root" or cleaned_taxonomy[0].lower() == "rootrank" or cleaned_taxonomy[0].lower() == "r:root":
        cleaned_taxonomy = cleaned_taxonomy[1:]
    # Complete taxonomy for uparse db
    if cleaned_taxonomy[0].startswith("d:"):
        tmp_tax = list()
        rank_idx = 0
        ranks = ["d:", "p:", "c:", "o:", "f:", "g:"]
        for taxa in cleaned_taxonomy:
            while not taxa.startswith(ranks[rank_idx]) and taxa != "Multi-affiliation" and taxa != "unclassified":
                tmp_tax.append(ranks[rank_idx] + "unknown_taxa")
                rank_idx += 1
            tmp_tax.append(taxa)
            rank_idx += 1
        cleaned_taxonomy = tmp_tax
    return cleaned_taxonomy

class SampleParameter(argparse.Action):
    """
    @summary : Argparse parameter for samples.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        # Set parser
        sample = getattr(namespace, self.dest)
        if sample is None:
            sample = dict()
        # Retrieve params
        for option in values:
            id, filepath = option.split(":")
            sample[id] = filepath
        setattr(namespace, self.dest, sample)


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Write a BIOM from grinder count profile.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-s', '--samples', type=str, action=SampleParameter, metavar=("SAMPLE_NAME:SAMPLE_PATH"), nargs='+', help="Samples names and grinder rank files." )
    group_input.add_argument( '-a', '--affiliation', required=True, help='Path to the databank source for simulated sequence (format: fasta). The description of sequences must be the taxonomy.' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-o', '--output', required=True, help='The output BIOM (format: BIOM).' )
    args = parser.parse_args()

    taxonomy_key = "real_taxonomy"
    biom = Biom( generated_by="grinder", matrix_type="sparse" )

    # Set observations count
    for sample_name in args.samples:
        biom.add_sample( sample_name )
        fh_abund = open( args.samples[sample_name] )
        for line in fh_abund: # Content format: "# rank<TAB>seq_id<TAB>rel_abund_perc"
            if not line.startswith('#'):
                fields = line.strip().split()
                try:
                    biom.add_observation( fields[1] )
                except: # already exist
                    pass
                biom.change_count( fields[1], sample_name, int(float(fields[2])*100000000000000) )################## depend de la precision grinder
        fh_abund.close()

    # Set taxonomy metadata
    fh_classif = FastaIO( args.affiliation )
    for record in fh_classif:
        try:
            metadata = biom.get_observation_metadata( record.id )
            if metadata is None or not metadata.has_key( taxonomy_key ):
                taxonomy = getCleanedTaxonomy(record.description)
                biom.add_metadata( record.id, taxonomy_key, taxonomy, "observation" )
        except ValueError: # is not in BIOM
            pass
    fh_classif.close()

    # Write BIOM
    BiomIO.write( args.output, biom )
