#!/usr/bin/env python2.7
#
# Copyright (C) 2019 INRA
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

__author__ = 'Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inra.fr'
__status__ = 'beta'

import re
import sys
import os
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
os.environ['PATH'] = CURRENT_DIR + os.pathsep + os.environ['PATH']

# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(CURRENT_DIR)), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import BiomIO


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
        ranks = ["d:", "p:", "c:", "o:", "f:", "g:","s:"]
        for taxa in cleaned_taxonomy:
            while not taxa.startswith(ranks[rank_idx]) and taxa != "Multi-affiliation" and taxa != "unclassified":
                tmp_tax.append(ranks[rank_idx] + "unknown_taxa")
                rank_idx += 1
            tmp_tax.append(taxa)
            rank_idx += 1
        cleaned_taxonomy = tmp_tax
    return cleaned_taxonomy

def get_realTax( taxonomy_key, input_biom ):
    """
    @summary: Returns count by taxa by rank in sample.
    @param input_biom: [str] Path to BIOM file.
    @return: [dict] The dictionary of count by taxa in dictionary by rank.
    """
    tax_list = list()
    biom = BiomIO.from_json( input_biom )
    
    for observation in biom.get_observations():
        taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
        if not taxonomy_clean in tax_list:
            tax_list.append(";".join(taxonomy_clean))
            
    return tax_list
    
def selectOneMultiaffiliation(real_tax, observation_id, possible_taxonomies):
    
    nb_select = 0
    selected = None
    for tax in possible_taxonomies:
        if tax in real_tax:
            selected = tax
            nb_select += 1
            
    if nb_select == 1:
        return selected
    else:
        if nb_select == 2:
            print "WARN : " + observation_id + " has multiple real correspondances"
        return possible_taxonomies[0]
    
def get_checkedTax( real_tax, input_biom, taxonomy_key, multi_affiliation ):
    """
    @summary:
    @param real_tax: [dict] Taxonomy.
    @param input_biom: [str] Path to BIOM file.
    @param taxonomy_key: [str] The metadata key for taxonomy.
    @param multi_affiliation: [bool] ************************************************************************************
    @return: [dict] The dictionary of count by taxa in dictionary by rank.   
    """
    tax_list = list()
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        # Get taxonomy
        if not multi_affiliation: # Standard affiliation
            taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
        else: # Multi-affiliation
            possible_taxonomies = [getCleanedTaxonomy(affi["taxonomy"]) for affi in observation["metadata"]["blast_affiliations"]]
            taxonomy_clean = selectOneMultiaffiliation(real_tax, observation["id"], possible_taxonomies)
        if taxonomy_clean not in tax_list:
            tax_list.append(";".join(taxonomy_clean))

    return tax_list

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Compare number of affiliations expected and detected.' )
    parser.add_argument( '-rk', '--real-tax-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in real biom. [Default: taxonomy]" )
    parser.add_argument( '-ck', '--checked-tax-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in checked biom. [Default: taxonomy]" )
    parser.add_argument( '-t', '--taxonomic-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels used in the metadata taxonomy.' )
    parser.add_argument( '-m', '--multi-affiliations', action='store_true', help='The taxonomy is produced by FROGS multi-affiliations.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-r', '--real-biom', required=True, help='Path to the theoretical abundance file (format: BIOM).' )
    group_input.add_argument( '-f', '--checked-biom', required=True, help='Path to the checked abundance file (format: BIOM).' )
    args = parser.parse_args()

    real_tax = get_realTax(args.real_tax_key, args.real_biom)
    checked_tax = get_checkedTax(real_tax, args.checked_biom, args.checked_tax_key, args.multi_affiliations)

    print "#Expected_tax\tDetected_tax\tRetrieved_tax"
    print str(len(real_tax)) + "\t" + str(len(checked_tax)) + "\t" + str(len(set(real_tax).intersection(checked_tax)))
    
