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
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import copy
import argparse
from collections import OrderedDict

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsBiom import *
from frogsSequenceIO import FastaIO

###################################################################################################################
###                                           FUNCTIONS                                                         ###
###################################################################################################################
def get_tax_consensus( taxonomies ):
    """
    @summary: Returns a consensus taxonomy from list of taxonomies.
    @param taxonomies: [list] The taxonomies to process. Each taxonomy is a list of rank taxon.
    @return: [list] The consensus taxonomy. The ambiguous ranks are replaced by "Multi-affiliation".
    @note:
        taxonomies = [ ["Bacteria", "Proteobacteria", "Gamma Proteobacteria", "Enterobacteriales"],
                       ["Bacteria", "Proteobacteria", "Beta Proteobacteria", "Methylophilales"] ]
        return = ["Bacteria", "Proteobacteria", "Multi-affiliation", "Multi-affiliation"]
    """

    consensus = list()
    if len(taxonomies) != 0:
        consensus = copy.copy(taxonomies[0])

    for curr_taxonomy in taxonomies[1:]:
        for rank, taxon in enumerate(curr_taxonomy):
            if consensus[rank] != "Multi-affiliation" and consensus[rank] != taxon:
                consensus[rank] = "Multi-affiliation"
    # Clean case with same taxon name in different branches:
    #      with taxonomies = [["A", "B", "C"], ["A", "L", "C"]]
    #      consensus is ["A", "Multi-affiliation", "C"] but must be ["A", "Multi-affiliation", "Multi-affiliation"]
    ancestor = ""
    for rank in range(len(consensus)):
        if ancestor == "Multi-affiliation":
            consensus[rank] = "Multi-affiliation"
        ancestor = consensus[rank]
    return consensus

def process(params):

    biom_in = BiomIO.from_json( params.input_biom )
    # check if biom_in has blast_taxonomy affiliations
    if not biom_in.has_metadata("blast_affiliations"):
        raise Exception("\n\n#ERROR : Your input biom file, "+ os.path.basename(params.input_biom) + ", does not contain any blast_affiliations metadata.\n\n")

    biom_out = Biom( generated_by='FROGS_aggregate_affiliated_otu', matrix_type="sparse" )

    # add samples in biom_out
    for sample_name in biom_in.get_samples_names():
        biom_out.add_sample( sample_name )

    # parse biom from most abondant OTU to less abondant one
    # save taxonomy
    # add OTU to biom_out if taxonomy is with poor %id or %cov or taxonomy not already saved
    # aggregate OTU to previous one if %id or %cov is big enough and share taxonomy with previous one

    # compute observation sum
    otu_sums = {}
    for otu_name,count_sum in biom_in.get_observations_counts():
        otu_sums[otu_name] = count_sum

    # save "confident" taxonomy
    otu_by_tax = dict()
    # save aggregated_otu_composition
    aggregated_otu = OrderedDict()
    otu_in = 0
    otu_out = 0
    otu_aggregated = 0

    # parse otu from most abondant to less ones
    for otu_name in sorted(otu_sums, key=lambda i: int(otu_sums[i]), reverse = True):
        otu_in += 1
        observation = biom_in.get_observations_by_name(otu_name)

        # is this OTU poorly affiliated
        min_id = 100
        min_cov = 100
        tax = list()
        for affiliation in observation["metadata"]["blast_affiliations"] : 
            if params.taxon_ignored and any(t in ";".join(affiliation["taxonomy"]) for t in params.taxon_ignored):
                continue
            if not affiliation["taxonomy"] in tax:
                tax.append(affiliation["taxonomy"])
            percent_id = affiliation["perc_identity"]
            percent_cov = affiliation["perc_query_coverage"]
            if percent_id < min_id : 
                min_id = percent_id
            if percent_cov < min_cov : 
                min_cov = percent_cov

        # Add otu because of poor affiliations stat
        if min_id < params.identity or min_cov < params.coverage :
            otu_out += 1
            biom_out.add_observation( otu_name, observation["metadata"] )
            for sample_name in biom_in.get_samples_names():
                count = biom_in.get_count(otu_name,sample_name)
                biom_out.add_count(otu_name, sample_name, count)
            aggregated_otu[otu_name] = list()
        # for confident taxonomy
        else:
            # check if all taxonomies are new
            is_new_tax = True
            equivalent_otu_name = ""

            for taxonomy in tax:
                if isinstance(taxonomy,list):
                    taxonomy = ";".join(taxonomy)
                if taxonomy in otu_by_tax:
                    is_new_tax = False
                    if equivalent_otu_name == "":
                        equivalent_otu_name = otu_by_tax[taxonomy]
                    elif otu_by_tax[taxonomy] != equivalent_otu_name:
                        Logger.static_write(params.log_file, '\tWarning: observation ' + otu_name + ' shares taxonomy ( '+ taxonomy +' with an other OTU : ' + otu_by_tax[taxonomy] + ', first detected OTU will be kept : ' + equivalent_otu_name + '\n' )

            # if new tax, add OTU and save taxonomies
            if is_new_tax:
                otu_out += 1
                biom_out.add_observation( otu_name, observation["metadata"] )
                for sample_name in biom_in.get_samples_names():
                    count = biom_in.get_count(otu_name, sample_name)
                    if count > 0 :
                        biom_out.add_count(otu_name, sample_name, count)
                aggregated_otu[otu_name] = list()
                for taxonomy in tax:
                    if isinstance(taxonomy,list):
                        taxonomy = ";".join(taxonomy)
                    otu_by_tax[taxonomy] = otu_name
            # else aggregation of OTU
            else:
                otu_aggregated += 1
                equivalent_otu = biom_out.get_observations_by_name(equivalent_otu_name)
                # add blast_affiliations
                aggregated_blast_affi = equivalent_otu["metadata"]["blast_affiliations"] + observation["metadata"]["blast_affiliations"]
                biom_out.add_metadata( equivalent_otu_name, "blast_affiliations", aggregated_blast_affi , subject_type="observation", erase_warning=False)
                # update consensus tax
                consensus_tax = get_tax_consensus([ affi["taxonomy"] for affi in aggregated_blast_affi ])
                biom_out.add_metadata( equivalent_otu_name, "blast_taxonomy", consensus_tax , subject_type="observation", erase_warning=False)
                # update counts
                for sample_name in biom_in.get_samples_names():
                    count = biom_out.get_count(equivalent_otu_name, sample_name) + biom_in.get_count(otu_name,sample_name)
                    biom_out.change_count(equivalent_otu_name, sample_name, count)
                # save aggregated composition
                aggregated_otu[equivalent_otu_name].append(otu_name)
                # update known taxonomies
                for taxonomy in tax:
                    if isinstance(taxonomy,list):
                        taxonomy = ";".join(taxonomy)
                    if not taxonomy in otu_by_tax:
                        otu_by_tax[taxonomy] = equivalent_otu_name

    # write biom output file
    BiomIO.write( params.output_biom, biom_out )

    # update fasta
    FH_in = FastaIO(params.input_fasta)
    FH_out = FastaIO(params.output_fasta, "wt")
    for record in FH_in:
        if record.id in aggregated_otu:
            FH_out.write(record)
    FH_in.close()
    FH_out.close()

    # write otu composition
    FH_compo = open(params.output_compo, "wt")
    for OTU in aggregated_otu:
        FH_compo.write(OTU + " " + " ".join(aggregated_otu[OTU]) + "\n")
    FH_compo.close()

    # simple log stat
    Logger.static_write(params.log_file, "# nb OTU in : "+str(otu_in) + "\n")
    Logger.static_write(params.log_file, "# nb OTU out : "+str(otu_out) + "\n")
    Logger.static_write(params.log_file, "# nb OTU aggregated : "+str(otu_aggregated) + "\n")

###################################################################################################################
###                                             MAIN                                                            ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Refine affiliations by aggregating OTU that share taxonomic affiliation with at least I% identity and C% coverage")
    parser.add_argument( '-i', '--identity', default=99.0, type=float, help="Min percentage identity to agggregate OTU. [Default: %(default)s]")
    parser.add_argument( '-c', '--coverage', default=99.0, type=float, help="Min percentage coverage to agggregate OTU. [Default: %(default)s]")
    parser.add_argument( '-t', '--taxon-ignored', type=str, nargs='*', help="Taxon list to ignore when OTUs agggregation")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-b', '--input-biom', required=True, help='Abundance table with affiliations metadata from the affiliation_OTU program (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help='OTU seed sequence file (format: Fasta).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', default='refined_affiliation.biom', help='File whith refind affiliation annotations. [Default: %(default)s]')
    group_output.add_argument('--output-compo', default='aggregated_otu_composition.tsv', help='Aggregated OTU composition [Default: %(default)s]')
    group_output.add_argument('--output-fasta', default='refined_affiliation.fasta', help='Updated OTU fasta file [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)
