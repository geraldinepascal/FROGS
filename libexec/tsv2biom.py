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

__author__ = 'Maria Bernard - Sigenae INRA'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'dev'

import os
import sys
import copy
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import Biom, BiomIO
from frogsSequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def store_multihits(input_multihits):
    multi_hit_dict=dict()
    FH_in = open(input_multihits)
    header_line = FH_in.readline()
    #to insure colnames transpositions in multiAffiFromBiom.py
    # from 1.1.0 :
    v1_1_0 = ["OTU", "Subject_taxonomy", "Blast_subject", "Prct_identity", "Prct_query_coverage", "e-value", "Alignment_length"]
    # to 1.1.2 (same as in TSV file):
    v1_2_0 = ["observation_name", "blast_taxonomy", "blast_subject", "blast_perc_identity", "blast_perc_query_coverage", "blast_evalue", "blast_aln_length"]
    for i in range(0,len(v1_1_0)):
        header_line = header_line.replace(v1_1_0[i],v1_2_0[i])

    header_line = header_line.replace("#","").strip().split()

    for colname in header_line:
        if colname not in ["observation_name", "blast_taxonomy", "blast_subject", "blast_perc_identity", "blast_perc_query_coverage", "blast_evalue", "blast_aln_length"] :
            raise Exception("\n"+colname+" is not a valid FROGS TSV column name.\nPlease restore the default column names: "+", ".join(["observation_name", "blast_taxonomy", " blast_subject", "blast_perc_identity", "blast_perc_query_coverage", "blast_evalue", "blast_aln_length"])+"\n\n")

    for line in FH_in.readlines():
        line = line.strip().split("\t")
        d=dict(zip(header_line,line))
        observation_name = d.pop("observation_name")
        if observation_name in multi_hit_dict:
            multi_hit_dict[observation_name].append(d)
        else:
            multi_hit_dict[observation_name] = [d]

    FH_in.close()

    return multi_hit_dict

def header_line_dict(fields,header,samples):
    meta_dict = dict()
    sample_dict = dict()
    seq_index=-1
    for idx,val in enumerate(header.strip().split()): 
        if val == "seed_sequence" :
            seq_index = idx
        elif val == "observation_sum":
            continue
        elif val in fields:
            meta_dict[idx]=val
        elif  val in samples:
            sample_dict[idx]=val

    return seq_index,meta_dict,sample_dict

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

def observation_blast_parts( metadata, obs_mutli_blast):
    """
    @summary: format blast_affiliation multi hit and compute consensus tax
    @param metadata: [dict] of metadata for one cluster keys = TSV fields.
    @param obs_multi_blast : [list] of dict for each blast hit for one cluster.
    @return blast_affiliations [list] of dict. Blast affiliations are filtered of non consistent taxonomy with metadata["blast_taxonomy"].
            consensus tax : computed from consitent blast_taxonomy and blast_hits.  
    """
    blast_taxonomy = ";".join([t for t in metadata["blast_taxonomy"] if t != "Multi-affiliation" ])
    blast_affiliations=[]

    # research of blast hit consistent with blast_taxonomy
    for hit in obs_mutli_blast:
        if blast_taxonomy in hit["blast_taxonomy"]:
            blast_affiliations.append({key.replace("blast_",""):hit[key] for key in hit})

    # compute consensus tax from blast_affiliations
    consensus = get_tax_consensus([ t["taxonomy"].split(";") for t in blast_affiliations])

    return consensus, blast_affiliations

def tsv_to_biom( input_tsv, multi_hit_dict, fields, samples_names, output_biom, output_fasta ):
    """
    @summary: Convert TSV file to Biom file.
    @param input_tsv: [str] Path to the TSV file.
    @param multi_hit_dict: [dict] Dictionnary describing equivalent multi blast hit : 
    dict[observation_name]=[ {"blast_taxonomy":taxonomy, "blast_subject":subject, "blast_perc_identity": per_id, "blast_perc_query_coverage":per_cov, "blast_evalue":eval, "blast_aln_length":aln}]
    @param fields: [list] column name to include as metadata (must at least contain observation_name): observation_sum and seed_sequence will be excluded, rdp_tax_and_bootstrap will be split in two metadata
    @param samples_names: [list] list of sample names.
    @param output_biom: [str] Path to the output file (format : BIOM).
    @param output_fasta: [str] Path to the output file (format : fasta).
    """
#     biom = Biom( generated_by='frogs', matrix_type="sparse" )
    biom = Biom( matrix_type="sparse" )

    seed_seq_idx = -1 
    metadata_index = dict()
    sample_index = dict()
    clusters_count = dict()
    clusters_metadata = dict()
    in_fh = open( input_tsv )

    if not output_fasta is None:
        Fasta_fh=FastaIO(output_fasta , "w" )

    # parse header and store column index 
    header=in_fh.readline()
    if header.startswith("#"):
        header=header[1:]
    header = header.strip()
    seed_seq_idx, metadata_index, sample_index = header_line_dict(fields,header,samples_names)
    if not output_fasta is None and seed_seq_idx == -1:
        raise Exception("\nYou want to extract seed fasta sequence but there is no seed_sequence column in your TSV file\n\n")

    # count by sample, and metadata
    for line in in_fh:

        cluster_name=""
        line_list=line.strip().split("\t")
        count_by_sample = {}
        metadata_dict = {}
        # parse columns
        for idx,val in enumerate(line_list):
            # recover metadata
            if idx in metadata_index:
                if metadata_index[idx]=="observation_name" :
                    cluster_name = val
                else:
                    metadata_dict[metadata_index[idx]] = val
            # recover samples count
            elif idx in sample_index and val > 0:
                count_by_sample[sample_index[idx]] = int(val)
            # recover seed sequence
            elif idx == seed_seq_idx:
                seed_seq = val

        # if fasta output file => store de seed sequence
        if not output_fasta is None:
            seq = Sequence( cluster_name, seed_seq) 
            Fasta_fh.write(seq)

        if "taxonomy" in metadata_dict:
            metadata_dict["taxonomy"] = metadata_dict["taxonomy"].split(";")

        # format rdp taxonomy to fit BIOM format
        if "rdp_tax_and_bootstrap" in metadata_dict:
            metadata_dict["rdp_taxonomy"]=[]
            metadata_dict["rdp_bootstrap"]=[]
            tax = metadata_dict["rdp_tax_and_bootstrap"].rstrip(";").split(";")
            for i in range(0,len(tax),2):
                metadata_dict["rdp_taxonomy"].append(tax[i])
                metadata_dict["rdp_bootstrap"].append(tax[i+1].replace("(","").replace(")",""))
            metadata_dict.pop("rdp_tax_and_bootstrap")

        # format blast taxonomy to fit BIOM format (one consensus blast_taxonomy and possible multiples blast_affiliation detailed
        if "blast_taxonomy" in metadata_dict:
            metadata_dict["blast_taxonomy"] = metadata_dict["blast_taxonomy"].split(";")

            # check multihit blast : filter non consistent taxonomy hit with blast_taxonomy (if TSV modified), and compute consensus tax (if multihit line suppressed)
            if metadata_dict["blast_subject"] == "multi-subject" and not multi_hit_dict is None:
                if not cluster_name in multi_hit_dict:
                    raise Exception("\n"+cluster_name+" has multi-subject tag but is not present in your multi-hit TSV file. Please, provide the original multi-hit TSV file.\n\n")
                else:
                    metadata_dict["blast_taxonomy"], metadata_dict["blast_affiliations"] = observation_blast_parts(metadata_dict, multi_hit_dict[cluster_name])
                    if metadata_dict["blast_affiliations"] == []:
                        raise Exception("\nyour multihit TSV file is no more consistent with your abundance TSV file for (at least) "+cluster_name+"\n\n")
            # no multi tag= blast affiliation is equal to blast_taxonomy
            else:
                blast_dict={key.replace("blast_",""):metadata_dict[key] for key in metadata_dict if key.startswith("blast")}
                metadata_dict["blast_affiliations"]=[blast_dict]

            # filter blast metadata which are moved to blast_affiliations
            for metadata in metadata_dict["blast_affiliations"][0]:
                if not metadata == "taxonomy":
                    metadata_dict.pop("blast_"+metadata)

        # add cluster and count to clusters_count dict
        clusters_count[cluster_name] = count_by_sample
        # ok print clusters_count[cluster_name].keys(), "CDT0#LOT05" in clusters_count[cluster_name], "CDT0#LOT02" in clusters_count[cluster_name]
        # add cluster and metadata to clusters_metadata dict
        clusters_metadata[cluster_name] = metadata_dict

    if not output_fasta is None:
        Fasta_fh.close()
    in_fh.close()

    #add samples to biom
    for sample_name in samples_names:
        biom.add_sample( sample_name )

    # add to cluster to biom
    for cluster_name in clusters_count:
        biom.add_observation( cluster_name, clusters_metadata[cluster_name] )
        for sample_name in samples_names:
            if clusters_count[cluster_name][sample_name] > 0:
                biom.add_count( cluster_name, sample_name, clusters_count[cluster_name][sample_name] )

    # Write
    BiomIO.write( output_biom, biom )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Convert TSV file to BIOM file.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--fields', default=['observation_name'], nargs='+', help="Metadata columns names (to include in the biom, you must at least have observation_name'. rdp_tax_and_bootstrap will be split in two taxonomy and bootstrap metadata, seed_sequence and observation_sum will be excluded.")
    parser.add_argument( '-s', '--samples-names', nargs='+', required=True, help="samples-names to include in the biom output")
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file', required=True, help='Path to the TSV file.' )
    group_input.add_argument( '-m', '--input-multihits', required=False, help='Path to the TSV multi hits file.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='Path to the output file (format : BIOM).')
    group_input.add_argument( '-a', '--output-fasta', default=None, required=False, help='Path to the output FASTA file.' )
    args = parser.parse_args()


    # Process
    multi_hit_dict = None
    if not args.input_multihits is None : 
        multi_hit_dict = store_multihits(args.input_multihits)
    tsv_to_biom( args.input_file, multi_hit_dict, args.fields, args.samples_names, args.output_file, args.output_fasta )
    