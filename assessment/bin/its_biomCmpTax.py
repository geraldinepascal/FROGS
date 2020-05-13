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

__author__ = 'Maria Bernard - SIGENAE Jouy en Josas'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
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
    first_rank = "d:"
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
    if cleaned_taxonomy[0].startswith("k:"):
        first_rank = "k:"
    if cleaned_taxonomy[0].startswith(first_rank):
        tmp_tax = list()
        rank_idx = 0
        ranks = [first_rank, "p:", "c:", "o:", "f:", "g:","s:"]
        for taxa in cleaned_taxonomy:
            while not taxa.startswith(ranks[rank_idx]) and taxa != "Multi-affiliation" and taxa != "unclassified" and taxa != "NA":
                tmp_tax.append(ranks[rank_idx] + "unknown_taxa")
                rank_idx += 1
            tmp_tax.append(taxa)
            rank_idx += 1
        while rank_idx != len(ranks):
            tmp_tax.append(ranks[rank_idx] + "unknown_taxa")
            rank_idx += 1
        
        cleaned_taxonomy = tmp_tax
    return cleaned_taxonomy

def getRealAbunByRank( taxonomy_key, input_biom, sample ):
    """
    @summary: Returns count by taxa by rank in sample.
    @param taxonomy_by_ref_id: [dict] Taxonomies by reference IDs.
    @param input_biom: [str] Path to BIOM file.
    @param sample: [str] sample name.
    @return: [dict] The dictionary of count by taxa in dictionary by rank.
    """
    abund_by_rank = list()
    tax_list = list()
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        count = biom.get_count( observation["id"], sample )
        if count > 0:
            taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
            if not ";".join(taxonomy_clean) in tax_list:
                tax_list.append(";".join(taxonomy_clean))
            for depth in range(len(taxonomy_clean)):
                if len(abund_by_rank) < depth+1:
                    abund_by_rank.append(dict())
                taxon = ";".join( taxonomy_clean[:depth+1] ) # prevent bug with same sp name but with different ancestors
                if not abund_by_rank[depth].has_key(taxon):
                    abund_by_rank[depth][taxon] = 0
                abund_by_rank[depth][taxon] += count
    return tax_list, abund_by_rank

def selectOneMultiaffiliation(real_tax, observation_id, possible_taxonomies, logfile):
    
    if len(possible_taxonomies) == 0:
        raise Exception ("\n" + observation_id + " has no blast_affiliations!\n")
        
    selected = list()
    for tax in possible_taxonomies:
        if ";".join(tax) in real_tax and ";".join(tax) not in selected:
            selected.append(";".join(tax))
            
    if len(selected) > 0:
        if len(selected) > 1:
            FH = open(logfile,"a")
            FH.write("WARN : " + observation_id + " has multiple real correspondances: \n")
            for tax in selected:
                FH.write("\t" + tax + "\n")
            FH.close()
        return selected[0].split(";")
    else :
        FH = open(logfile,"a")
        FH.write("WARN : " + observation_id + " has no real correspondances, return first multiaffiliation\n")
        FH.close()
        return possible_taxonomies[0]
    
    
def getCheckedAbunByRank( real_tax, input_biom, sample, taxonomy_key, multi_affiliation, logfile):
    """
    @summary:
    @param real_tax: [dict] Taxonomy by reference IDs.
    @param input_biom: [str] Path to BIOM file.
    @param sample: [str] sample name.
    @param taxonomy_key: [str] The metadata key for taxonomy.
    @param multi_affiliation: [bool] ************************************************************************************
    @return: [dict] The dictionary of count by taxa in dictionary by rank.   
    """
    abund_by_rank = list()
    tax_list = list()
    full_tax_list = list()
    nb_seq = 0
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        count = biom.get_count( observation["id"], sample )
        if count > 0:
            nb_seq += 1 
            # Get taxonomy
            if not multi_affiliation: # Standard affiliation
                taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
                if ";".join(taxonomy_clean) not in full_tax_list:
                        full_tax_list.append(";".join(taxonomy_clean))
            else: # Multi-affiliation
                possible_taxonomies = [getCleanedTaxonomy(affi["taxonomy"]) for affi in observation["metadata"]["blast_affiliations"]]
                for taxonomy_clean in possible_taxonomies:
                    if ";".join(taxonomy_clean) not in full_tax_list:
                        full_tax_list.append(";".join(taxonomy_clean))
                taxonomy_clean = selectOneMultiaffiliation(real_tax, observation["id"], possible_taxonomies, logfile)
            
            if ";".join(taxonomy_clean) not in tax_list:
                tax_list.append(";".join(taxonomy_clean))

            # Store count
            for depth in range(len(taxonomy_clean)):
                if len(abund_by_rank) < depth+1:
                    abund_by_rank.append(dict())
                taxon = ";".join( taxonomy_clean[:depth+1] ) # prevent bug with same sp name but with different ancestors
                if not abund_by_rank[depth].has_key(taxon):
                    abund_by_rank[depth][taxon] = 0
                abund_by_rank[depth][taxon] += count
    return nb_seq, full_tax_list, tax_list, abund_by_rank


def cmpTaxAbund( expected, obtained, depth ):
    """
    @warning: same depth for the 2 files
    """
    identity = 0
    divergence = 0
    common_taxa = 0
    expected_specific = 0
    obtained_specific = 0
    already_processed = list()
    total_expected = sum([abund for taxon, abund in expected[depth].items()])
    total_obtained = sum([abund for taxon, abund in obtained[depth].items()])
    for taxon, abund in expected[depth].items():
        already_processed.append( taxon )
        if obtained[depth].has_key(taxon):
            common_taxa += 1
            prct_expected = (float(abund)*100)/total_expected
            prct_obtained = (float(obtained[depth][taxon])*100)/total_obtained
            identity += min( prct_expected, prct_obtained )
            if abs(prct_expected - prct_obtained) > 1:
                print "DIFF\t" + str(abs(prct_expected - prct_obtained)) + "\t" + taxon + "\t" + str(prct_expected) + "\t" + str(prct_obtained)
        else:
            expected_specific += 1
            if (float(abund)*100)/total_expected > 1:
                print "GRINDER\t" + str((float(abund)*100)/total_expected) + "\t" + taxon
    for taxon, abund in obtained[depth].items():
        if not taxon in already_processed:
            obtained_specific += 1
            if (float(abund)*100)/total_expected > 1:
                print "CHECKED\t" + str((float(abund)*100)/total_expected) + "\t" + taxon
    divergence = 100 - identity
    return {
        'divergence': divergence,
        'common_taxa': common_taxa,
        'expected_specific': expected_specific,
        'obtained_specific': obtained_specific,
    }


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Compare expected affiliations to produced affiliations for selected sample.' )
    parser.add_argument( '-s', '--sample', type=str, required=True, help="The name of the selected sample." )
    parser.add_argument( '-rk', '--real-tax-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in real biom. [Default: taxonomy]" )
    parser.add_argument( '-ck', '--checked-tax-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in checked biom. [Default: taxonomy]" )
    parser.add_argument( '-t', '--taxonomic-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels used in the metadata taxonomy.' )
    parser.add_argument( '-m', '--multi-affiliations', action='store_true', help='The taxonomy is produced by FROGS multi-affiliations.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-r', '--real-biom', required=True, help='Path to the theoretical abundance file (format: BIOM).' )
    group_input.add_argument( '-f', '--checked-biom', required=True, help='Path to the checked abundance file (format: BIOM).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-l', '--log-file', required=True, help='Path to log file.' )
    args = parser.parse_args()

    #decimal_precision = "%.5f"

    real_tax , real_tax_abundance = getRealAbunByRank( args.real_tax_key, args.real_biom, args.sample )

    nb_seq , all_checked_tax, checked_tax, checked_tax_abundance = getCheckedAbunByRank( real_tax, args.checked_biom, args.sample, args.checked_tax_key, args.multi_affiliations, args.log_file )
    
    print "#Expected_tax\tCluster\tDetected_tax\tRetrieved_tax"
    print str(len(real_tax)) + "\t" + str(nb_seq) + "\t" + str(len(all_checked_tax)) + "\t" + str(len(set(real_tax).intersection(all_checked_tax)))
    
    print ""

    print "#Rank\tDivergence (%)\tCommon\tReal specific\tChecked specific"
    for depth in range(len(checked_tax_abundance)):
        metrics = cmpTaxAbund( real_tax_abundance, checked_tax_abundance, depth )
        print str(args.taxonomic_ranks[depth]) + "\t" + str(metrics['divergence']) + "\t" + str(metrics['common_taxa']) + "\t" + str(metrics['expected_specific']) + "\t" + str(metrics['obtained_specific'])
