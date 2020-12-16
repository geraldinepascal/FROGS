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
__version__ = '0.2.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'beta'

import re
import sys
import json
import argparse
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
        ranks = ["d:", "p:", "c:", "o:", "f:", "g:"]
        for taxa in cleaned_taxonomy:
            while not taxa.startswith(ranks[rank_idx]) and taxa != "Multi-affiliation" and taxa != "unclassified":
                tmp_tax.append(ranks[rank_idx] + "unknown_taxa")
                rank_idx += 1
            tmp_tax.append(taxa)
            rank_idx += 1
        cleaned_taxonomy = tmp_tax
    return cleaned_taxonomy


def getDuplicationGroupByID( input_file ):
    """
    @summary: 
    @param input_file: [str] Path to the file with group of duplicated sequences by line.
    @returns: [dict] By duplicated sequence ID the IDs list of identical sequence.
    """
    group_by_id = dict()
    fh_input = open(input_file)
    for line in fh_input:
        group = line.strip().split("\t")
        for group_member in group:
            group_by_id[group_member] = group
    fh_input.close()
    return group_by_id


def getRealTaxByRefID( input_biom, taxonomy_key, duplication_groups ):
    """
    @summary: Return taxonomy by reference.
    @param input_biom: [str] Path to BIOM file.
    @param taxonomy_key: [str] The metadata key for taxonomy.
    @param duplication_groups: [dict] By reference ID the list of references with the same sequence.
    @return: [dict] List of taxonomies by reference ID.
             Example: 
               {
                 "MVF01000012.1.1317": [
                   ["Root", "Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacteriales", "Enterobacteriaceae", "Cronobacter", "Escherichia coli BIDMC 73"]
                 ],
                 "JQ607252.1.1437": [
                   ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus", "bacterium NLAE-zl-P471"],
                   ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus", "Staphylococcus aureus M17299"]
                 ] 
               }
    """
    taxonomy_by_obs_id = dict()
    tmp_taxonomy_by_obs_id = dict()
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
        taxonomy_by_obs_id[observation["id"]] = [taxonomy_clean]
        tmp_taxonomy_by_obs_id[observation["id"]] = taxonomy_clean
    if duplication_groups is not None:
        for obs_id in duplication_groups:
            taxonomy_by_obs_id[obs_id] = list()
            for id_duplicated_seq in duplication_groups[obs_id]: # For each duplication group member
                taxonomy_by_obs_id[obs_id].append(tmp_taxonomy_by_obs_id[id_duplicated_seq])
    return taxonomy_by_obs_id


def getRealAbunByRank( taxonomy_by_ref_id, input_biom, sample ):
    """
    @summary: Returns count by taxa by rank in sample.
    @param taxonomy_by_ref_id: [dict] Taxonomies by reference IDs.
    @param input_biom: [str] Path to BIOM file.
    @param sample: [str] sample name.
    @return: [dict] The dictionary of count by taxa in dictionary by rank.
    """
    abund_by_rank = list()
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        count = biom.get_count( observation["id"], sample )
        if count > 0:
            taxonomy_clean = taxonomy_by_ref_id[observation["id"]][0]
            for depth in range(len(taxonomy_clean)):
                if len(abund_by_rank) < depth+1:
                    abund_by_rank.append(dict())
                taxon = ";".join( taxonomy_clean[:depth+1] ) # prevent bug with same sp name but with different ancestors
                if not abund_by_rank[depth].has_key(taxon):
                    abund_by_rank[depth][taxon] = 0
                abund_by_rank[depth][taxon] += count
    return abund_by_rank

def refIDIsRetrieved( ref_id, subjects_ids, duplication_groups ):
    ref_id_is_retrieved = False
    if duplication_groups is None or ref_id not in duplication_groups: # Sequence is uniq in dataset
        if ref_id in subjects_ids:
            ref_id_is_retrieved = True
    else: # Sequence has duplicated sequence in dataset
        for id_duplicated_seq in duplication_groups[ref_id]:
            if id_duplicated_seq in subjects_ids:
                ref_id_is_retrieved = True
    return ref_id_is_retrieved


def taxIsRetrieved( ref_taxonomies, subjects_taxonomies ):
    tax_is_retrieved = False
    if len(ref_taxonomies) == 1: # Sequence is uniq in dataset
        if ";".join(ref_taxonomies[0]) in subjects_taxonomies:
            tax_is_retrieved = True
    else: # Sequence has duplicated sequence in dataset
        for current_ref_tax in ref_taxonomies:
            if ";".join(current_ref_tax) in subjects_taxonomies:
                tax_is_retrieved = True
    return tax_is_retrieved


def getCheckedAbunByRank( real_tax, input_biom, sample, taxonomy_key, multi_affiliation, duplication_groups ):
    """
    @summary:
    @param real_tax: [dict] Taxonomy by reference IDs.
    @param input_biom: [str] Path to BIOM file.
    @param sample: [str] sample name.
    @param taxonomy_key: [str] The metadata key for taxonomy.
    @param multi_affiliation: [bool] ************************************************************************************
    @param duplication_groups: [dict] By reference ID the list of IDs for references with the same sequence.
    @return: [dict] The dictionary of count by taxa in dictionary by rank.   
    """
    abund_by_rank = list()
    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        count = biom.get_count( observation["id"], sample )
        if count > 0:
            # Get taxonomy
            ref_id = observation["metadata"]["grinder_source"]
            taxonomy_clean = getCleanedTaxonomy(observation["metadata"][taxonomy_key])
            if not multi_affiliation: # Standard affiliation
                if not "," in ref_id: # Non chimera
                    if taxIsRetrieved(real_tax[ref_id], [taxonomy_clean]):
                        taxonomy_clean = real_tax[ref_id][0]
            else: # Multi-affiliation
                if not "," in ref_id: # Non chimera
                    subjects_ids = [affi["subject"] for affi in observation["metadata"]["blast_affiliations"]]
                    possible_taxonomies = [";".join(getCleanedTaxonomy(affi["taxonomy"])) for affi in observation["metadata"]["blast_affiliations"]]
                    # Manage ambiguity
                    if refIDIsRetrieved(ref_id, subjects_ids, duplication_groups):
                        taxonomy_clean = real_tax[ref_id][0]
                    elif len(subjects_ids) > 499 and taxIsRetrieved(real_tax[ref_id], possible_taxonomies):
                        taxonomy_clean = real_tax[ref_id][0]
                    elif "Multi-affiliation" in taxonomy_clean:
                        taxonomy_clean = getCleanedTaxonomy(observation["metadata"]["blast_affiliations"][0]["taxonomy"]) # Select one
                else: # Chimera
                    if "Multi-affiliation" in taxonomy_clean:
                        taxonomy_clean = getCleanedTaxonomy(observation["metadata"]["blast_affiliations"][0]["taxonomy"]) # Select one
            # Store count
            for depth in range(len(taxonomy_clean)):
                if len(abund_by_rank) < depth+1:
                    abund_by_rank.append(dict())
                taxon = ";".join( taxonomy_clean[:depth+1] ) # prevent bug with same sp name but with different ancestors
                if not abund_by_rank[depth].has_key(taxon):
                    abund_by_rank[depth][taxon] = 0
                abund_by_rank[depth][taxon] += count
    return abund_by_rank


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
    group_input.add_argument( '-u', '--uniq-groups', default=None, help='Path to the file with by line the list of IDs of initial sequences with the same sequence (format: TSV).' )
    args = parser.parse_args()

    #decimal_precision = "%.5f"
    duplic_groups_by_id = None
    if args.uniq_groups is not None:
        duplic_groups_by_id = getDuplicationGroupByID(args.uniq_groups)
    real_taxonomy_by_ref_id = getRealTaxByRefID( args.real_biom, args.real_tax_key, duplic_groups_by_id )
    real_tax_abundance = getRealAbunByRank( real_taxonomy_by_ref_id, args.real_biom, args.sample )
    checked_tax_abundance = getCheckedAbunByRank( real_taxonomy_by_ref_id, args.checked_biom, args.sample, args.checked_tax_key, args.multi_affiliations, duplic_groups_by_id )

    print "#Rank\tDivergence (%)\tCommon\tReal specific\tChecked specific"
    for depth in range(len(checked_tax_abundance)):
        metrics = cmpTaxAbund( real_tax_abundance, checked_tax_abundance, depth )
        print str(args.taxonomic_ranks[depth]) + "\t" + str(metrics['divergence']) + "\t" + str(metrics['common_taxa']) + "\t" + str(metrics['expected_specific']) + "\t" + str(metrics['obtained_specific'])
