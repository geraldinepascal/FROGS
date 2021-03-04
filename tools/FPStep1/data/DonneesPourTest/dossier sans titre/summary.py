#!/usr/bin/env python2.7
# -*-coding:Utf-8 -*

#Import

import os
import sys
import argparse
import json
import re
from numpy import median
#from cmd import Cmd
    """
    @summary: Writes the process summary in one html file.
    @param summary_file: [str] path to the output html file.
    @param file_fasta: [str] path to the fasta file of OTU
    @param excluded_tree: [str] path to the Newick file.
    @param file_tree: [str] path to the Newick file.

    """

 # to summary FASTA number && abundances number               
    summary_info = {
       'sequence_kept' : 0,
       'sequence_removed' : 0     
    }
    number_seq_all = 0
    # to detail removed OTU
    removed_details_categories =["Taxonomic Information", "Sequence Number", "% with abundance total", "Sequence length"]
    removed_details_data =[]
    
    # to build one metadata for tree view
    dic_seq={}
    list_seq_all=list()
    list_out_tree=[]
    """
    biom=BiomIO.from_json(biomfile)
    treefile = open(treefile, "r")
    newick = treefile.read().strip()
    """

    # record nb OTU and abundance
    for otu in FastaIO(fasta_in):
        list_seq_all.append(seq.id)
        number_seq_all +=1
        #number_abundance_all += biom.get_observation_count(otu.id)



    # record details about removed OTU
    if align_out is not None:
        for otu in FastaIO(align_out):
            summary_info['otu_removed'] +=1
            summary_info['abundance_removed'] += biom.get_observation_count(otu.id)
            
            # to built one table of OTUs out of phylogenetic tree
            taxonomy=""
            if biom.has_metadata("taxonomy"):
                taxonomy = ";".join(biom.get_observation_metadata(otu.id)["taxonomy"]) if issubclass(biom.get_observation_metadata(otu.id)["taxonomy"].__class__,list) else str(biom.get_observation_metadata(otu.id)["taxonomy"])
            elif biom.has_metadata("blast_taxonomy"): 
                taxonomy = ";".join(biom.get_observation_metadata(otu.id)["blast_taxonomy"]) if issubclass(biom.get_observation_metadata(otu.id)["blast_taxonomy"].__class__,list) else str(biom.get_observation_metadata(otu.id)["blast_taxonomy"])
            abundance=biom.get_observation_count(otu.id)
            percent_abundance=abundance*100/(float(number_abundance_all))
            length=len(otu.string)
            info={"name": otu.id, "data": [taxonomy, abundance, percent_abundance, length]}
            removed_details_data.append(info)
            list_out_tree.append(otu.id)
