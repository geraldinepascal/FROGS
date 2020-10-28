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
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse
import copy

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsSequenceIO import *
from frogsUtils import *
from frogsBiom import *

###################################################################################################################
###                                           FUNCTIONS                                                         ###
###################################################################################################################
def get_fasta_size(fasta_file):
    """
    @summary : return dictionnary of sequence lengh : dict[ref_id] = size
    """
    size_by_ref = {}    
    FH_in = FastaIO(fasta_file)
    for record in FH_in:
        size_by_ref[record.id] = len(record.string)

    return size_by_ref

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

def select_smallest(blast_affiliations, ref_size):
    """
    @summary : reduce the list of blast affiliations be keeping the smallest reference
    @param blast_affiliations : [list] list of dictionnary 
                    [{
                        "perc_identity": 100.0, 
                        "taxonomy": [
                            "Bacteria", 
                            "Bacteroidetes", 
                            "Flavobacteriia", 
                            "Flavobacteriales", 
                            "Flavobacteriaceae", 
                            "Pibocella", 
                            "Pibocella ponti"
                        ], 
                        "evalue": "0.0", 
                        "aln_length": 421, 
                        "perc_query_coverage": 100.0, 
                        "subject": "AY576654.1.1447"
                    }]
    """
    min_size = 0
    smallest_affi = []
    all_affi = [affi["subject"] for affi in blast_affiliations]

    for affi in blast_affiliations:
        if not affi["subject"] in ref_size:
            return blast_affiliations
        if min_size == 0 or ref_size[affi["subject"]] < min_size :
            smallest_affi = [affi]
            min_size = ref_size[affi["subject"]]
        elif ref_size[affi["subject"]] == min_size:
            smallest_affi.append(affi)

    return smallest_affi

def process(params) :
    """
    @summary : Select smallest amplicon reference among multiaffiliations
    """

    # save reference fasta file
    ref_size = get_fasta_size(params.ITS_reference)

    # biom initialization
    biom_in = BiomIO.from_json( params.input_biom )
    biom_out = Biom( matrix_type="sparse" )

    if not biom_in.has_metadata("blast_taxonomy"):
        raise Exception("\n\n#ERROR : Your biom file need to be affiliated with FROGS_affiliation_OTU tool\n\n")

    #add samples to biom
    for sample_name in biom_in.get_samples_names():
        biom_out.add_sample( sample_name )

    nb_obs = 0
    nb_obs_multi_affi = 0
    nb_obs_multi_affi_resolved = 0
    nb_multi_affi_removed = 0
    for observation in biom_in.get_observations():
        nb_obs += 1
        # reduce multiaffiliations list
        if len(observation['metadata']["blast_affiliations"]) > 1:
            nb_obs_multi_affi +=1
            new_blast_affi = select_smallest(observation['metadata']['blast_affiliations'], ref_size)
            if len(new_blast_affi) < len(observation['metadata']['blast_affiliations']):
                nb_multi_affi_removed += (len(observation['metadata']["blast_affiliations"]) - len(new_blast_affi))
                nb_obs_multi_affi_resolved += 1
            observation['metadata']['blast_affiliations'] = new_blast_affi
            consensus_tax = get_tax_consensus( [ affi["taxonomy"] for affi in new_blast_affi])
            observation['metadata']['blast_taxonomy'] = consensus_tax

        # add observation in biom_out
        biom_out.add_observation( observation['id'], observation['metadata'])
        for sample_name in biom_in.get_samples_names():
            if biom_in.get_count(observation['id'], sample_name) > 0:
                biom_out.add_count( observation['id'], sample_name, biom_in.get_count(observation['id'], sample_name))    

    # simple log stat
    Logger.static_write(params.log_file, "# nb OTU : "+str(nb_obs) + "\n")
    Logger.static_write(params.log_file, "# nb OTU with multiaffiliations : "+str(nb_obs_multi_affi) + "\n")
    Logger.static_write(params.log_file, "# nb OTU with multiaffiliations resolved : "+str(nb_obs_multi_affi_resolved) + "\n")
    Logger.static_write(params.log_file, "# nb multiaffiliations removed : "+str(nb_multi_affi_removed) + "\n")

    BiomIO.write( params.output_biom, biom_out )
###################################################################################################################
###                                           MAIN                                                              ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Refine affiliations, to manage ITS amplicon included in other ITS sequence.")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-b', '--input-biom', required=True, help='Abundance table with affiliations metadata from the affiliation_OTU program (format: BIOM).')
    group_input.add_argument('-r', '--ITS-reference', required=True, help='reference ITS1 or ITS2 fasta file')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-biom', default='refined_affiliation.biom', help='File whith refind affiliation annotations. [Default: %(default)s]')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)
