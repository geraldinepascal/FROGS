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
        raise_exception( Exception("\n\n#ERROR : Your input biom file, "+ os.path.basename(params.input_biom) + ", does not contain any blast_affiliations metadata.\n\n"))

    biom_out = Biom( generated_by='FROGS_aggregate_affiliated_asv', matrix_type="sparse" )

    # add samples in biom_out
    for sample_name in biom_in.get_samples_names():
        biom_out.add_sample( sample_name )

    # parse biom from most abondant ASV to less abondant one
    # save taxonomy
    # add ASV to biom_out if taxonomy is with poor %id or %cov or taxonomy not already saved
    # aggregate ASV to previous one if %id or %cov is big enough and share taxonomy with previous one

    # compute observation sum
    asv_sums = {}
    for asv_name,count_sum in biom_in.get_observations_counts():
        asv_sums[asv_name] = count_sum

    # save "confident" taxonomy
    asv_by_tax = dict()
    # save aggregated_asv_composition
    aggregated_asv = OrderedDict()
    asv_in = 0
    asv_out = 0
    asv_aggregated = 0

    # parse asv from most abondant to less ones
    for asv_name in sorted(asv_sums, key=lambda i: int(asv_sums[i]), reverse = True):
        asv_in += 1
        observation = biom_in.get_observations_by_name(asv_name)

        # is this ASV poorly affiliated
        min_id = 100
        min_cov = 100
        tax = list()
        for affiliation in observation["metadata"]["blast_affiliations"] : 
            if params.taxon_ignored and any(t in ";".join(affiliation["taxonomy"]) for t in params.taxon_ignored):
                continue
            if affiliation["perc_identity"] == "no data" or affiliation["perc_query_coverage"] == "no data": # Check if it is not an unaffiliated ASV
                continue
            if not affiliation["taxonomy"] in tax:
                tax.append(affiliation["taxonomy"])
            percent_id = float(affiliation["perc_identity"])
            percent_cov = float(affiliation["perc_query_coverage"])
            if percent_id < min_id : 
                min_id = percent_id
            if percent_cov < min_cov : 
                min_cov = percent_cov

        # Add asv because of poor affiliations stat
        if min_id < params.identity or min_cov < params.coverage :
            asv_out += 1
            biom_out.add_observation( asv_name, observation["metadata"] )
            for sample_name in biom_in.get_samples_names():
                count = biom_in.get_count(asv_name,sample_name)
                biom_out.add_count(asv_name, sample_name, count)
            aggregated_asv[asv_name] = list()
        # for confident taxonomy
        else:
            # check if all taxonomies are new
            is_new_tax = True
            equivalent_asv_name = ""

            for taxonomy in tax:
                if isinstance(taxonomy,list):
                    taxonomy = ";".join(taxonomy)
                if taxonomy in asv_by_tax:
                    is_new_tax = False
                    if equivalent_asv_name == "":
                        equivalent_asv_name = asv_by_tax[taxonomy]
                    elif asv_by_tax[taxonomy] != equivalent_asv_name:
                        Logger.static_write(params.log_file, '\tWarning: observation ' + asv_name + ' shares taxonomy ( '+ taxonomy +' with an other ASV : ' + asv_by_tax[taxonomy] + ', first detected ASV will be kept : ' + equivalent_asv_name + '\n' )

            # if new tax, add ASV and save taxonomies
            if is_new_tax:
                asv_out += 1
                biom_out.add_observation( asv_name, observation["metadata"] )
                for sample_name in biom_in.get_samples_names():
                    count = biom_in.get_count(asv_name, sample_name)
                    if count > 0 :
                        biom_out.add_count(asv_name, sample_name, count)
                aggregated_asv[asv_name] = list()
                for taxonomy in tax:
                    if isinstance(taxonomy,list):
                        taxonomy = ";".join(taxonomy)
                    asv_by_tax[taxonomy] = asv_name
            # else aggregation of ASV
            else:
                asv_aggregated += 1
                equivalent_asv = biom_out.get_observations_by_name(equivalent_asv_name)
                # add blast_affiliations
                aggregated_blast_affi = equivalent_asv["metadata"]["blast_affiliations"] + observation["metadata"]["blast_affiliations"]
                biom_out.add_metadata( equivalent_asv_name, "blast_affiliations", aggregated_blast_affi , subject_type="observation", erase_warning=False)
                # update consensus tax
                consensus_tax = get_tax_consensus([ affi["taxonomy"] for affi in aggregated_blast_affi ])
                biom_out.add_metadata( equivalent_asv_name, "blast_taxonomy", consensus_tax , subject_type="observation", erase_warning=False)
                # update counts
                for sample_name in biom_in.get_samples_names():
                    count = biom_out.get_count(equivalent_asv_name, sample_name) + biom_in.get_count(asv_name,sample_name)
                    biom_out.change_count(equivalent_asv_name, sample_name, count)
                # save aggregated composition
                aggregated_asv[equivalent_asv_name].append(asv_name)
                # update known taxonomies
                for taxonomy in tax:
                    if isinstance(taxonomy,list):
                        taxonomy = ";".join(taxonomy)
                    if not taxonomy in asv_by_tax:
                        asv_by_tax[taxonomy] = equivalent_asv_name

    # write biom output file
    BiomIO.write( params.output_biom, biom_out )

    # update fasta
    FH_in = FastaIO(params.input_fasta)
    FH_out = FastaIO(params.output_fasta, "wt")
    for record in FH_in:
        if record.id in aggregated_asv:
            FH_out.write(record)
    FH_in.close()
    FH_out.close()

    # write asv composition
    FH_compo = open(params.output_compo, "wt")
    for ASV in aggregated_asv:
        FH_compo.write(ASV + " " + " ".join(aggregated_asv[ASV]) + "\n")
    FH_compo.close()

    # simple log stat
    Logger.static_write(params.log_file, "# nb ASV in : "+str(asv_in) + "\n")
    Logger.static_write(params.log_file, "# nb ASV out : "+str(asv_out) + "\n")
    Logger.static_write(params.log_file, "# nb ASV aggregated : "+str(asv_aggregated) + "\n")

###################################################################################################################
###                                             MAIN                                                            ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Refine affiliations by aggregating ASV that share taxonomic affiliation with at least I% identity and C% coverage")
    parser.add_argument( '-i', '--identity', default=99.0, type=float, help="Min percentage identity to agggregate ASV. [Default: %(default)s]")
    parser.add_argument( '-c', '--coverage', default=99.0, type=float, help="Min percentage coverage to agggregate ASV. [Default: %(default)s]")
    parser.add_argument( '-t', '--taxon-ignored', type=str, nargs='*', help="Taxon list to ignore when ASVs agggregation")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-b', '--input-biom', required=True, help='Abundance table with affiliations metadata from the taxonomic_affiliation tool (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help='ASV seed sequence file (format: Fasta).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', default='refined_affiliation.biom', help='File whith refind affiliation annotations. [Default: %(default)s]')
    group_output.add_argument('--output-compo', default='aggregated_asv_composition.tsv', help='Aggregated ASV composition [Default: %(default)s]')
    group_output.add_argument('--output-fasta', default='refined_affiliation.fasta', help='Updated ASV fasta file [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)
