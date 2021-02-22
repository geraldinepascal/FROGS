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

__author__ = 'Plateforme bioinformatique Toulouse / Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'


import re
import sys
import argparse
import warnings
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
        
    # Deal with OTUs without affiliation
    if len(cleaned_taxonomy) == 0:
		cleaned_taxonomy = ["d:unknown_taxa", "p:unknown_taxa", "c:unknown_taxa", "o:unknown_taxa", "f:unknown_taxa", "g:unknown_taxa","s:unknown_taxa"]
        
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


def cmpTaxAbund( expected, checked, depth ):
    identity = 0
    divergence = 0
    common_taxa = 0
    expected_specific = 0
    checked_specific = 0
    detailed_cmp = dict()
    already_processed = list()
    total_expected = sum([abund for taxon, abund in expected[depth].items()])
    
    #~ if len(checked) == 0:
        #~ a=1
        #~ return {
        #~ 'divergence': 100.0,
        #~ 'common_taxa': 0,
        #~ 'expected_specific': 0,
        #~ 'checked_specific': checked_specific,
        #~ 'detailed_cmp': detailed_cmp
    #~ }
    #~ else:
		#~ min_level_in_checked = max(checked.keys()) # 6 for Species, 5 for Genus (case where no Species were detected)

    total_checked = sum([abund for taxon, abund in checked[depth].items()])
    for taxon, abund in expected[depth].items():
        already_processed.append( taxon )
        if checked[depth].has_key(taxon):
            common_taxa += 1
            prct_expected = (float(abund)*100)/total_expected
            prct_checked = (float(checked[depth][taxon])*100)/total_checked
            identity += min( prct_expected, prct_checked )
            detailed_cmp[taxon] = {"expected": prct_expected, "checked": prct_checked}
        else:
            expected_specific += 1
            prct_expected = (float(abund)*100)/total_expected
            detailed_cmp[taxon] = {"expected": prct_expected, "checked": 0}
    for taxon, abund in checked[depth].items():
        if not taxon in already_processed:
            checked_specific += 1
            prct_checked = (float(checked[depth][taxon])*100)/total_checked
            detailed_cmp[taxon] = {"expected": 0, "checked": prct_checked}
    divergence = 100 - identity
    return {
        'divergence': divergence,
        'common_taxa': common_taxa,
        'expected_specific': expected_specific,
        'checked_specific': checked_specific,
        'detailed_cmp': detailed_cmp
    }


def get_expected( abund_file ):
    expected_by_depth = dict()
    FH_expected = open(abund_file)
    for line in FH_expected:
        taxonomy, count = line.strip().split("\t")
        clean_taxonomy = getCleanedTaxonomy(taxonomy)
        for rank_depth in range(len(clean_taxonomy)):
            rank_taxonomy = ";".join(clean_taxonomy[:rank_depth + 1])
            if rank_depth not in expected_by_depth:
                expected_by_depth[rank_depth] = dict()
            if rank_taxonomy not in expected_by_depth[rank_depth]:
                expected_by_depth[rank_depth][rank_taxonomy] = 0
            expected_by_depth[rank_depth][rank_taxonomy] += float(count)
    FH_expected.close()
    return expected_by_depth


def get_checked( abund_file, checked_sample, taxonomy_key, expected_by_depth ):
    
    checked_by_depth = dict()
    biom = BiomIO.from_json(abund_file)
    
    for current_obs in biom.get_observations():
        
        #recuperation de la taxo
        clean_taxonomy = getCleanedTaxonomy(current_obs["metadata"][taxonomy_key]) if current_obs["metadata"][taxonomy_key] is not None else ["unknown_taxa"]*len(expected_by_depth)
        #print(clean_taxonomy)

        # recuperation du count
        count = biom.get_count(current_obs["id"], checked_sample)
        if count > 0:
            # check multiaffi
            if clean_taxonomy[len(clean_taxonomy)-1] == "Multi-affiliation":
                #print(">>>> ",clean_taxonomy)
                selected = list()
                taxonomies = list()
                expected_taxonomies = expected_by_depth[len(clean_taxonomy)-1]
                # pour chaque multi-affi
                for affi_idx in range(len(current_obs["metadata"]["blast_affiliations"])):
                    # clean taxo au format string
                    affi_taxonomy = ";".join(getCleanedTaxonomy(current_obs["metadata"]["blast_affiliations"][affi_idx]["taxonomy"]))
                    #print("#####################", affi_taxonomy)
                    # compte des taxo qui sont attendues
                    if affi_taxonomy not in taxonomies:
                        taxonomies.append(affi_taxonomy)
                        if affi_taxonomy in expected_taxonomies:
                            selected.append((getCleanedTaxonomy(current_obs["metadata"]["blast_affiliations"][affi_idx]["taxonomy"])))
                            #print(affi_taxonomy," ????????????????????")
                #print("\t", selected)
                if len(selected) > 0:
					#clean_taxonomy = selected[0]
					for multi in selected:
						for rank_depth in range(len(multi)):
							rank_taxonomy = ";".join(multi[:rank_depth + 1])
							if rank_depth not in checked_by_depth:
								checked_by_depth[rank_depth] = dict()
							if rank_taxonomy not in checked_by_depth[rank_depth]:
								checked_by_depth[rank_depth][rank_taxonomy] = 0
							checked_by_depth[rank_depth][rank_taxonomy] += count/len(selected)
                    
                    #if len(selected) > 1 :
                        #warnings.warn( "Multi-affiliation cannot be resolved for " + str((float(count)*100)/biom.get_total_count()) + "% sequences. Possible taxonomies: '" + "', '".join(taxonomies) + "'." )
                        
                else:
                    clean_taxonomy = taxonomies[0].split(";")
                    for rank_depth in range(len(clean_taxonomy)):
						rank_taxonomy = ";".join(clean_taxonomy[:rank_depth + 1])
						if rank_depth not in checked_by_depth:
							checked_by_depth[rank_depth] = dict()
						if rank_taxonomy not in checked_by_depth[rank_depth]:
							checked_by_depth[rank_depth][rank_taxonomy] = 0
						checked_by_depth[rank_depth][rank_taxonomy] += count
                    
                
                #return checked_by_depth
            #print("\n")
            #print(clean_taxonomy)
            #print("\n")
            else:
				for rank_depth in range(len(clean_taxonomy)):
					rank_taxonomy = ";".join(clean_taxonomy[:rank_depth + 1])
					if rank_depth not in checked_by_depth:
						checked_by_depth[rank_depth] = dict()
					if rank_taxonomy not in checked_by_depth[rank_depth]:
						checked_by_depth[rank_depth][rank_taxonomy] = 0
					checked_by_depth[rank_depth][rank_taxonomy] += count
    
    return checked_by_depth


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Compare expected abundance by taxon (in INPUT_TSV) with obtained abudance for the sample (sample in INPUT_BIOM).' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-e', '--expected-abund', required=True, help='Path to the expected abundance by taxon (format: TSV).' )
    parser.add_argument( '-c', '--checked-abund', required=True, help='Path to the checked abundance file (format: BIOM).' )
    parser.add_argument( '-s', '--checked-sample', required=True, help='Name of checked sample.' )
    parser.add_argument( '-k', '--taxonomy-key', default="taxonomy", help='The used tag to store taxonomy in the BIOM file. [Default: taxonomy].' )
    parser.add_argument( '-r', '--checked-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='Evaluated ranks.' )
    args = parser.parse_args()

    # Expected
    expected_by_depth = get_expected(args.expected_abund)
    
    # Retrieved
    checked_by_depth = get_checked(args.checked_abund, args.checked_sample, args.taxonomy_key, expected_by_depth)
    #print checked_by_depth
    #for i in sorted(expected_by_depth[6]):
	#	print i, expected_by_depth[6][i]
    #print "\n------\n"
    #for i in sorted(checked_by_depth[6]):
	#	print i, checked_by_depth[6][i]
    #sys.exit()
    # Comparison
    details = dict()
    print "#Rank\tDivergence (%)\tCommon\tExpected specific\tChecked specific"
    #try:
    for depth, rank in enumerate(args.checked_ranks):
        metrics = cmpTaxAbund( expected_by_depth, checked_by_depth, depth )
        #print metrics
        details[depth] = metrics["detailed_cmp"]
        print rank + "\t" + str(metrics['divergence']) + "\t" + str(metrics['common_taxa']) + "\t" + str(metrics['expected_specific']) + "\t" + str(metrics['checked_specific'])
    print ""
    #except KeyError:
    #print rank + "\t" + str(metrics['divergence']) + "\t" + str(metrics['common_taxa']) + "\t" + str(metrics['expected_specific']) + "\t" + str(metrics['checked_specific'])
        #print rank + "\t" + str(100.0) + "\t" + str(0) + "\t" + str(len(expected_by_depth[len(args.checked_ranks)-1])) + "\t" + str(0)

    try:
        for depth, rank in enumerate(args.checked_ranks):
            print "#Rank " + rank
            print "#Taxon\tExpected (%)\tChecked (%)"
            for taxon in sorted(details[depth]):
                print taxon + "\t" + str(details[depth][taxon]["expected"]) + "\t" + str(details[depth][taxon]["checked"])
            print ""
    except KeyError:
        pass
