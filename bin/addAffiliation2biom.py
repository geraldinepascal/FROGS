#!/usr/bin/env python2.7
#
# Copyright (C) 2014 INRA
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
__version__ = '2.2.1'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import sys
import copy
import argparse
from biom import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

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

def get_tax_from_fasta( reference_file ):
    """
    @summary: Returns taxonomies by sequence ID from a reference databank.
    @param reference_file: [str] Path to the fasta with taxonomy in sequences descriptions.
    @return: [dict] The taxonomies by sequence ID.
    @note: Fasta sequence header example : '>ADC62191 Bacteria; Proteobacteria; Gammaproteobacteria; Chromatiales; Chromatiaceae; Allochromatium.; Allochromatium vinosum DSM 180'
    """
    taxonomy_by_seqid = {}
    for line in open(reference_file).readlines():
        if line.startswith('>'):
            parts = line[1:].strip().split(None, 1)  # Line example : '>ADC62191 Bacteria; Proteobacteria; Gammaproteobacteria; Chromatiales; Chromatiaceae; Allochromatium.; Allochromatium vinosum DSM 180'
            if parts[1].endswith(';'):
                parts[1] = parts[1][:-1] # Removes last ';'
            taxonomy = list()
            for rank, taxon in enumerate(parts[1].split(";")):
                taxon_name = taxon.split('[id:')[0].strip() # Removes '[id:XXX]' in taxon name
                if rank != 0 or taxon_name != "Root": # Skips the taxon named 'Root'.
                    taxonomy.append(taxon_name)
            taxonomy_by_seqid[parts[0]] = taxonomy
    return taxonomy_by_seqid

def get_rdp_affi( rdp_files ):
    """
    @summary: Returns taxonomies and bootstraps by query ID from a RDP result files.
    @param rdp_files: [list] Path to the RDP result files.
    @return: [dict] The taxonomies and bootstrap by sequence ID.
    """
    rdp_annot = dict()
    for current_rdp in rdp_files:
        FH_rdp = open(current_rdp)
        for line in FH_rdp:
            fields = line.strip().split("\t")
            cluster_id = fields[0]
            taxonomy = list()
            bootstraps = list()
            for idx in range(2, len(fields[2:])):
                if (2 - idx) % 3 == 0 and fields[idx + 1] != "rootrank":
                    taxonomy.append(fields[idx].split('[id:')[0].strip())
                    bootstraps.append(float(fields[idx + 2]))
            rdp_annot[cluster_id] = dict()
            rdp_annot[cluster_id]['taxonomy'] = taxonomy
            rdp_annot[cluster_id]['bootstrap'] = bootstraps
        FH_rdp.close()
    return rdp_annot

def get_bests_blast_affi( blast_files, taxonomy_by_subject ):
    """
    @summary: Returns by query ID the blast metrics for the best HSP(s).
    @param blast_files: [list] Path to the blast 'outfmt 6' result files.
    @param taxonomy_by_subject: [dict] The taxonomies by subject sequence ID.
    @return: [dict] The blast metrics by sequence ID. Example:
        [
            'query_12': [
                {
                    "perc_identity": 100.0,
                    "taxonomy": [
                        "Bacteria",
                        "Bacteroidetes",
                        "Bacteroidia",
                        "Bacteroidales",
                        "Bacteroidaceae",
                        "Bacteroides",
                        "unknown species"
                    ],
                    "evalue": "0.0",
                    "perc_query_coverage": 100.0,
                    "subject": "FJ507190.1.1382",
                    "aln_length": 450
                },
                {
                    "perc_identity": 100.0,
                    "taxonomy": [
                    "Bacteria",
                        "Bacteroidetes",
                        "Bacteroidia",
                        "Bacteroidales",
                        "Bacteroidaceae",
                        "Bacteroides",
                        "unknown species"
                    ],
                    "evalue": "0.0",
                    "perc_query_coverage": 100.0,
                    "subject": "FJ507188.1.1382",
                    "aln_length": 450
                },
            ],
            'query_15': [
            ...
        ]

    """
    blast_annot = dict()
    for current_blast in blast_files:
        FH_blast = open(current_blast)
        for line in FH_blast:
            parts = line.strip().split()
            query_id = parts[0].split(";")[0]
            score = float(parts[11])
            if not blast_annot.has_key(query_id) or blast_annot[query_id]['score'] < score:
                blast_annot[query_id] = {
                    'score': score,
                    'alignments': list(),
                }
            if blast_annot[query_id]['score'] == score: # select best HSP
                subject_id = parts[1].split("#")[0]  # Subject field : <ID>#<PARTIAL_DESC>
                blast_annot[query_id]['alignments'].append({
                    'subject': subject_id,
                    'taxonomy': taxonomy_by_subject[subject_id],
                    'evalue': parts[10],
                    'aln_length': int(parts[3]),
                    'perc_identity': float(parts[2]),
                    'perc_query_coverage': (int(parts[7]) - int(parts[6]) + 1) / float(parts[12]) * 100
                })
        FH_blast.close()
    return blast_annot

def aff_to_metadata(reference_file, biom_in, biom_out, blast_files=None, rdp_files=None):
    """
    @summary: Add taxonomy metadata on biom file from a blast result.
    @param reference_file: [str] The path to the reference file.
    @param biom_in: [str] The path to the Biom file to process.
    @param biom_out: [str] The path to the biom output file.
    @param blast_files: [list] the list of the path to the blast results in tabular format (outfmt 6 with NCBI Blast+).
    @param rdp_files: [list] the list of path to the RDPClassifier results.
    """
    # Build an hash with the taxonomy for each gene (key=gene_id ; value=gene_taxonomy)
    taxonomy_by_reference = get_tax_from_fasta( reference_file )

    # Retrieve blast clusters annotations
    cluster_blast_annot = dict()
    if blast_files is not None:
        cluster_blast_annot = get_bests_blast_affi( blast_files, taxonomy_by_reference )
    del taxonomy_by_reference

    # Retrieve rdp clusters annotations
    cluster_rdp_annot = dict()
    if rdp_files is not None:
        cluster_rdp_annot = get_rdp_affi( rdp_files )

    # Add metadata to biom
    biom = BiomIO.from_json(biom_in)
    for cluster in biom.get_observations():
        cluster_id = cluster["id"]
        # Blast
        if blast_files is not None:
            blast_taxonomy = None
            blast_affiliations = list()
            if cluster_blast_annot.has_key(cluster_id): # Current observation has a match
                blast_taxonomy = get_tax_consensus( [alignment['taxonomy'] for alignment in cluster_blast_annot[cluster_id]['alignments']] )
                blast_affiliations = cluster_blast_annot[cluster_id]['alignments']
            biom.add_metadata( cluster_id, "blast_affiliations", blast_affiliations, "observation" )
            biom.add_metadata( cluster_id, "blast_taxonomy", blast_taxonomy, "observation" )
        # RDP
        if rdp_files is not None:
            rdp_taxonomy = None
            rdp_bootstrap = None
            if cluster_rdp_annot.has_key(cluster_id):
                rdp_taxonomy = cluster_rdp_annot[cluster_id]['taxonomy']
                bootstrap = cluster_rdp_annot[cluster_id]['bootstrap']
            biom.add_metadata(cluster_id, "rdp_taxonomy", rdp_taxonomy, "observation")
            biom.add_metadata(cluster_id, "rdp_bootstrap", bootstrap, "observation")
    BiomIO.write(biom_out, biom)


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Add affiliation results from Blast and/or RDP to BIOM file.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )

    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument('-f', '--reference', required=True, help='Preformated reference file.')
    group_input.add_argument( '-b', '--blast-file', nargs="*", default=None, required=False, help='Path to Blast output file.' )
    group_input.add_argument( '-r', '--rdp-file', nargs="*", default=None, required=False, help='Path to RDP Classifier output file.' )
    group_input.add_argument( '-i', '--biom-in', default=None, required=True, help='Path to biom input file' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--biom-out', default=None, required=True, help='Path to biom output file' )
    args = parser.parse_args()

    # Process
    if args.blast_file is None and args.rdp_file is None:
        raise Exception("at least one blast or one RDPClassifier output file is needed\n")
    else:
        aff_to_metadata(args.reference, args.biom_in, args.biom_out, args.blast_file, args.rdp_file)