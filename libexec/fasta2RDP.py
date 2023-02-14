#!/usr/bin/env python3
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
#-*-coding:utf-8-*-
__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'dev'


import sys,os
import random
import argparse
import numpy as np
import unicodedata
from frogsNode import *
from frogsSequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_taxonomy( node ):
    """
    @summary: Returns the taxonomy of the node.
    @param node: [Node] The node processed.
    @return: [str] The taxonomy based on 'clean_name'. Each level is separated by ';'.
    """
    parent_taxonomy = ";".join([getCleanName(ancestor) for ancestor in node.get_ancestors()[1:]])
    if parent_taxonomy == "":
        return getCleanName(node)
    return parent_taxonomy + ";" + getCleanName(node)

def get_cleaned_sp( species_name ):
    """
    @summary: Returns the species name without strain information.
    @param species_name: [str] The standard name of the species (example: last field in taxonomy extract from the header fasta).
    @return: [str] The cleaned species name.
    """
    new_species_name = species_name
    if species_name.startswith("uncultured") or species_name.lower() in ["undefined", "unidentified", "incertae sedis"]:
        new_species_name = "unknown species"
    else:
        # Remove strain information
        pattern = re.compile( '^([^\s]+ sp\.).+' )
        matches = pattern.match( new_species_name )
        if matches is not None:
            new_species_name = matches.group(1)
        else:
            pattern = re.compile( '^([^\s]+ [^\s]+) subsp\.' )
            matches = pattern.match( new_species_name )
            if matches is not None:
                new_species_name = matches.group(1)
            else:
                pattern = re.compile( '(.+) DSMY? \d+$' )
                matches = pattern.match( new_species_name )
                if matches is not None:
                    new_species_name = matches.group(1)
    return new_species_name


def getCleanName(node):
    clean_name = node.name
    if node.name.lower() != "root":
        clean_name = node.name + " [id: " + str(node.metadata["id"]) + "]"
    return clean_name


def writeRDPTax( FH_tax, node ):
    """
    @summary: Writes the node and all the descendant. The Output format is RDP trainning database taxonomy (RDPTools/classifier.jar train).
    @param FH_tax: [File] The file handle on output file.
    @param node: [Node] The node to write.
    """
    node_depth = node.get_depth()
    parent_id = None
    if node_depth == 0:
        parent_id = -1
    else:
        parent_id = node.parent.metadata["id"]
    #taxid*taxon name*parent taxid*depth*rank
    FH_tax.write( str(node.metadata["id"]) + '*' + getCleanName(node) + '*' + str(parent_id) + '*' + str(node_depth) + '*' + str(node.metadata["rank"]) + "\n" )
    for child in node.get_children():
        writeRDPTax( FH_tax, child )


def treeFromFasta( in_fasta, ranks ):
    tree = Node("Root", None, None, {"id":0, "rank":"rootrank"})
    current_id = 1
    FH_databank = FastaIO(in_fasta)
    for record in FH_databank:
        desc = unicodedata.normalize('NFD', record.description).encode('ascii', 'ignore').decode('utf-8')
        if desc.endswith(";"):
            desc=desc[:-1]
        taxonomy = [taxa.strip() for taxa in desc.split(";")]
        parent = tree
        for rank_depth, taxa in enumerate(taxonomy):
            if not ranks is None:
                rank = ranks[rank_depth]
            else :
                rank=taxa[0]
            if not parent.has_child(taxa): 
            #####################################################################" pb niv espece car diff sp in connu du meme genre donne meme nom
                taxa_node = Node(taxa, parent, None, {"id":current_id, "rank":rank})
                current_id += 1
            parent = parent.get_child(taxa)
    FH_databank.close()
    return tree


def mergeFastaTax( in_tax, in_fasta, out_fasta ):
    tax_by_id = dict()
    # Get taxonomies by sequence ID
    FH_tax = open(in_tax)
    for line in FH_tax:
        seq_id, taxonomy = line.strip().split(None, 1)
        tax = unicodedata.normalize('NFD', taxonomy).encode('ascii', 'ignore').decode('utf-8')
        if tax.endswith(";"):
            tax_by_id[seq_id] = tax[:-1]
        else:
            tax_by_id[seq_id] = tax
    FH_tax.close()
    # Write fasta with taxonomy
    FH_in_db = FastaIO(in_fasta)
    FH_out_db = FastaIO(out_fasta, "wt")
    for record in FH_in_db:
        record.description = tax_by_id[record.id]
        FH_out_db.write(record)
    FH_in_db.close()
    FH_out_db.close()


def writeRDPFasta(tree, in_fasta, out_fasta):
    FH_in_db = FastaIO(in_fasta)
    FH_out_db = FastaIO(out_fasta, "wt")
    for record in FH_in_db:
        desc = unicodedata.normalize('NFD', record.description).encode('ascii', 'ignore').decode('utf-8')
        if desc.endswith(";"):
            desc=desc[:-1]
        taxonomy = [taxa.strip() for taxa in desc.split(";")]
        clean_taxonomy = [getCleanName(tree)]
        node = tree
        for taxa in taxonomy:
            tax = unicodedata.normalize('NFD', taxa).encode('ascii', 'ignore').decode('utf-8')
            node = node.get_child(tax)
            clean_taxonomy.append(getCleanName(node))
        record.description = ";".join(clean_taxonomy)
        FH_out_db.write(record)
    FH_in_db.close()
    FH_out_db.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Prepare file for RDP classifier trainings. (NB : input file(s) must be utf-8 encoding)")
    )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-r', '--ranks' , nargs='*', help='The ordered ranks levels used in the metadata taxonomy. [Default: first letter of each taxon (usefull for Greengenes)]' ) 
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-d', '--databank', required=True, help='The reference databank (format: FASTA). Each sequence must have the same number of tacxonomy level and the header must have this format: "ID<TAB>TAX_LVL1;TAX_LVL2". or provide the taxonomy with --taxonomy option' )
    group_input.add_argument( '-t', '--taxonomy', required=False, help='The reference databank (format: TSV). Each sequence must have the same number of taxonomy level and the header must have this format: "ID<TAB>TAX_LVL1;TAX_LVL2".' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--rdp-taxonomy', required=True, help='The selected sequences (format: RDPTax).' )
    group_output.add_argument( '--rdp-fasta', required=True, help='The selected sequences (format: RDPTax).' )
    args = parser.parse_args()

    # Pre-process
    db_with_tax = args.databank
    if args.taxonomy is not None:
        db_with_tax = args.databank+".tmp" #################################################
        mergeFastaTax(args.taxonomy, args.databank, db_with_tax)
    
    # Build tree
    databank_tree = treeFromFasta(db_with_tax, args.ranks)

    # Write RDP tax
    FH_RDPTax = open(args.rdp_taxonomy, "wt")
    writeRDPTax(FH_RDPTax, databank_tree)
    FH_RDPTax.close()

    # Write RDP fasta
    writeRDPFasta(databank_tree, db_with_tax, args.rdp_fasta)

	# remove tmp file
    if args.taxonomy is not None:
        os.remove(args.databank+".tmp")

