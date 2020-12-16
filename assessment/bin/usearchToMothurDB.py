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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import re
import argparse
from frogsSequenceIO import *


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


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Convert the databank provided by usearch in databank usable by mothur.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    parser.add_argument( '-d', '--databank', required=True, help='Path to the usearch databank (format: fasta). The ID of sequences must be the ended by taxonomy (format: "<ID>;tax=<TAX>").' )
    parser.add_argument( '-t', '--taxonomy', required=True, help='Path to the outputed taxonomy file (format: mothur tax).' )
    parser.add_argument( '-f', '--fasta', required=True, help='Path to the outputed sequences file (format: fasta). IDs are separated of the taxonomy.' )
    args = parser.parse_args()

    # Process
    FH_tax_out = open(args.taxonomy, "w")
    FH_fasta_out = FastaIO(args.fasta, "w")
    FH_databank_in = FastaIO(args.databank)
    for record in FH_databank_in:
        if record.description is None:
            if ";tax=" in record.id:
                record.id, record.description = record.id.split(";tax=")
            else:
                raise Exception("Taxonomy cannot be retrieved for '" + record.id + "'.")
        taxonomy = getCleanedTaxonomy(record.description)
        FH_tax_out.write(record.id + "\t" + ";".join(taxonomy) + ";\n" )
        record.description = ";".join(taxonomy) + ";"
        FH_fasta_out.write(record)
    FH_tax_out.close()
