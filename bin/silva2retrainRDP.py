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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.3.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'dev'


import re
import argparse
from node import Node


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
    parent_taxonomy = ";".join([ancestor.metadata["clean_name"] for ancestor in node.get_ancestors()[1:]])
    if parent_taxonomy == "":
        return node.metadata["clean_name"]
    return parent_taxonomy + ";" + node.metadata["clean_name"]

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

def write_rdp_tax( FH_tax, node ):
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
    FH_tax.write( str(node.metadata["id"]) + '*' + node.metadata["clean_name"] + '*' + str(parent_id) + '*' + str(node_depth) + '*' + node.metadata["rank"] + "\n" )
    for child in node.get_children():
        write_rdp_tax( FH_tax, child )

def silva_tax_2_tree( taxonomy_file, authorized_ranks, domains_filtered=None ):
    """
    @summary: Returns a taxonomic tree from the silva taxonomy.
    @param taxonomy_file: [str] The silva taxonomy file (provide taxon name, taxon level and taxon rank name). Line example: 'Archaea;Crenarchaeota;Thermoprotei;    7    class        119'
    @param authorized_ranks: [list] The ranks reported in tree.
    @param domains_filtered: [list] These domains are not kept in tree.
    @return: [Node, dict, int] The root node of the tree, the link between complete taxonomy (all levels) and node (only authorized ranks), the id of the next new element added in tree.
    """
    next_new_id = 1000000000
    taxonomy_ref = dict()
    taxonomy_tree = Node( "Root", None, None, {'clean_name':'Root', 'id':0, 'rank':'rootrank'} )
    FH_tax = open( taxonomy_file )
    for line in FH_tax:
        if line.strip() != "":
            # Parse line
            matches = re.search( "^(.+);\s+(\d+)\s+([^\s]+)", line.strip() )
            if matches is not None:
                evaluated_taxonomy = matches.group(1).split(";")
                evaluated_id = matches.group(2)
                evaluated_rank = matches.group(3)
            else:
                raise Exception("Incorrect line content : '" + line.strip() + "'.")
            if domains_filtered is None or not evaluated_taxonomy[0].strip().lower() in domains_filtered:
                # Go to the most downstream already existing ancestor of the element
                parent_node = taxonomy_tree
                if len(evaluated_taxonomy) > 1:
                    parent_node = taxonomy_ref[";".join(evaluated_taxonomy[:-1]).lower()]
                current_node = None
                # Add evaluated node if it has a valid rank
                if evaluated_rank in authorized_ranks:
                    if authorized_ranks.index(evaluated_rank) <= authorized_ranks.index(parent_node.metadata["rank"]):
                        if evaluated_taxonomy[-1] == "uncultured" and parent_node.name == evaluated_taxonomy[-2]:
                            evaluated_rank = authorized_ranks[authorized_ranks.index(parent_node.metadata["rank"])+1]
                        else:
                            raise Exception( "The taxonomy in file '" + taxonomy_file + "' seems to be incoherent. The taxon '" + ";".join(evaluated_taxonomy) + "' is tagged as '" + evaluated_rank + "' and its ancestor '" + get_taxonomy(parent_node) + "' is tagged as '" + parent_node.metadata["rank"] + "'." )
                    # Complete missing ranks between parent and evaluated
                    evaluated_rank_depth = authorized_ranks.index(evaluated_rank)
                    while authorized_ranks[evaluated_rank_depth -1] != parent_node.metadata["rank"]:
                        missing_rank_depth = authorized_ranks.index(parent_node.metadata["rank"]) +1
                        missing_name = "unknown " + authorized_ranks[missing_rank_depth]
                        missing_node = None
                        if parent_node.has_child( missing_name ):
                            missing_node = parent_node.get_child( missing_name )
                        else:
                            missing_node = Node( missing_name, parent_node, {}, {"clean_name" : missing_name + " [id:" + str(next_new_id) + "]", "id" : next_new_id, "rank" : authorized_ranks[missing_rank_depth]} )
                            next_new_id += 1
                            parent_node.add_child( missing_node )
                        parent_node = missing_node
                    # Add evaluated node
                    evaluated_name = evaluated_taxonomy[-1].strip()
                    if evaluated_name.lower() in ["unidentified", "uncultured", "incertae sedis", "unknown " + evaluated_rank, parent_node.name + " incertae sedis"]: # Clean unknown taxon name
                        evaluated_name = "unknown " + evaluated_rank
                    if evaluated_name == "unknown " + evaluated_rank and parent_node.has_child(evaluated_name): # Is unknown and unknown already exists
                        current_node = parent_node.get_child( evaluated_name )
                    else: # Is not unknown or is unknown but does not already exist
                        current_node = Node( evaluated_name, parent_node, {}, {"clean_name" : evaluated_name + " [id:" + evaluated_id + "]", "id" : evaluated_id, "rank" : evaluated_rank} )
                        parent_node.add_child( current_node )
                # Store link between complete taxonomy and node in tree
                if current_node is None:
                    current_node = parent_node
                taxonomy_ref[";".join(evaluated_taxonomy).lower()] = current_node
    FH_tax.close()

    return taxonomy_ref, taxonomy_tree, next_new_id

def process( input_fasta, input_tax, output_fasta, output_tax, domains_filtered=None, ranks=["rootrank", "domain", "phylum", "class", "order", "family", "genus", "species"] ):
    """
    @param input_fasta: [str] The silva sequences file.
    @param input_tax: [str] The silva taxonomy file.
    @param output_fasta: [str] The RDP retrain classifier sequence file.
    @param output_tax: [str] The RDP retrain classifier taxonomy file.
    @param domains_filtered: [list] These domains are not kept in output.
    @param ranks: [list] The ranks reported in outputs.
    """
    # Create tax tree
    taxonomy_ref, taxonomy_tree, new_taxon_id = silva_tax_2_tree( input_tax, ranks, domains_filtered )

    # Write fasta
    FH_cleanFasta = open( output_fasta, "w" )
    FH_fasta = open( input_fasta )
    is_filtered = None
    for line in FH_fasta:
        if line.startswith( ">" ):
            line_fields = line.strip().split()
            evaluated_id = line_fields[0]
            evaluated_taxonomy = " ".join(line_fields[1:]).split(";")
            if domains_filtered is not None and evaluated_taxonomy[0].strip().lower() in domains_filtered:
                is_filtered = True
            elif not taxonomy_ref.has_key( ";".join(evaluated_taxonomy[:-1]).lower() ):
                is_filtered = True
                print "The sequence '" + evaluated_id + "' is skipped because the node for '" + ";".join(evaluated_taxonomy) + "' does not exist in the taxonomy file."
            else:
                is_filtered = False
                clean_taxonomy = get_taxonomy(taxonomy_ref[";".join(evaluated_taxonomy[:-1]).lower()])
                if not "species" in ranks:
                    raise Exception( "The execution without 'species' rank is not implemented." )
                else:
                    parent_node = taxonomy_ref[";".join(evaluated_taxonomy[:-1]).lower()]
                    # Go to genus
                    while parent_node.metadata["rank"] != "genus":
                        parent_depth = ranks.index( parent_node.metadata["rank"] )
                        if parent_node.has_child( "unknown " + ranks[parent_depth +1] ):
                            parent_node = parent_node.get_child( "unknown " + ranks[parent_depth +1] )
                        else:
                            missing_name = "unknown " + ranks[parent_depth +1]
                            missing_node = Node( missing_name, parent_node, {}, {"clean_name" : missing_name + " [id:" + str(new_taxon_id) + "]", "id" : new_taxon_id, "rank" : ranks[parent_depth +1]} )
                            new_taxon_id += 1
                            parent_node.add_child( missing_node )
                            parent_node = missing_node
                    # Add species to tree
                    species_name = get_cleaned_sp(evaluated_taxonomy[-1])
                    species_node = None
                    if parent_node.has_child( species_name ):
                        species_node = parent_node.get_child( species_name )
                    else:
                        species_node = Node( species_name, parent_node, {}, {"clean_name": species_name.replace('*', ' ').replace('<', ' ').replace('>', ' ') + " [id:" + str(new_taxon_id) + "]", "id": new_taxon_id, "rank": "species" } )
                        new_taxon_id += 1
                        parent_node.add_child( species_node )
                    clean_taxonomy = get_taxonomy(species_node)
                FH_cleanFasta.write( line_fields[0] + '\tRoot' + ';' + clean_taxonomy + "\n" )
        elif not is_filtered:
            FH_cleanFasta.write( line.strip().replace('u', 't').replace('U', 'T') + "\n" )
    FH_fasta.close()
    FH_cleanFasta.close()

    # Write RDP tax file
    FH_cleanTax = open( output_tax, "w" )
    write_rdp_tax( FH_cleanTax, taxonomy_tree )
    FH_cleanTax.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=('''
Converts silva databank to RDP databank.
Command sequence to generate an databank usable with classifier.jar :
    python silva2retrainRDP.py --reference silva.fasta --taxonomy silva.tax --rdp-reference rdp_silva.fasta --rdp-taxonomy rdp_silva.tax
    java -Xmx60g -jar <RDP_TOOLS_PATH>/classifier.jar train -o . -s rdp_silva.fasta -t rdp_silva.tax
    cp <RDP_TOOLS_PATH>/classifier/samplefiles/rRNAClassifier.properties rdp_silva.fasta.properties''')
    )
    parser.add_argument( '-f', '--filtered', default=None, nargs='+', choices=['bacteria', 'archaea', 'eukaryota'], help='The filtered domains.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-r', '--reference', required=True, help='The silva sequences file (format : FASTA). Example : SILVA_119_SSURef_Nr99_tax_silva.fasta.' )
    group_input.add_argument( '-t', '--taxonomy', required=True, help='The silva taxonomy file (format : TSV). Example : tax_slv_ssu_nr_119.txt.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--rdp-reference', default='rdp_seq.fasta', help='The RDP sequence file (FASTA).' )
    group_output.add_argument( '--rdp-taxonomy', default='rdp_tax.txt', help='The RDP taxonomy file.' )
    args = parser.parse_args()

    # Process
    process( args.reference, args.taxonomy, args.rdp_reference, args.rdp_taxonomy, args.filtered )

"""
##################################################################################################################################################
#
# TEST DATASETS
#
##################################################################################################################################################
Dataset tree
***********
((((((("sp6")"genus4")"unknown family")"order4")"class4")"unknown phylum")"domain2",(((((("sp2")"unknown genus")"family2")"unknown order")"class2",(((("sp1")"genus1")"family1")"order1")"class1")"unknown phylum",((((("sp3")"genus2")"family3")"unknown order")"class3",(((("sp4")"genus3")"family4")"order3",((("sp5")"unknown genus")"unknown family")"unknown order")"unknown class")"phylum1")"domain1")"Root"

Dataset 1
*********
    - databank.fasta:
>A1 domain1;class1;order1;family1;genus1;sp1
A
>B1 domain1;class2;unidentified;family2;sp2
A
>C1 domain1;phylum1;class3;family3;genus2;sp3
A
>D1 domain1;phylum1;order3;family4;genus3;sp4
A
>E1 domain1;phylum1;sp5
A
>F1 domain2;class4;order4;genus4;sp6
A

    - databank.tax:
domain1;    1    domain    119
domain1;class1;    11    class    119
domain1;class1;order1;    12    order    119
domain1;class1;order1;family1;    13    family    119
domain1;class1;order1;family1;genus1;    14    genus    119
domain1;class2;    23    class    119
domain1;class2;unidentified;    24    order    119
domain1;class2;unidentified;family2;    25    family    119
domain1;phylum1;    26    phylum    119
domain1;phylum1;class3;    27    class    119
domain1;phylum1;class3;family3;    28    family    119
domain1;phylum1;class3;family3;genus2;    29    genus    119
domain1;phylum1;order3;    38    order    119
domain1;phylum1;order3;family4;    39    family    119
domain1;phylum1;order3;family4;genus3;    40    genus    119
domain2;    52    domain    119
domain2;class4;    53    class    119
domain2;class4;order4;    61    order    119
domain2;class4;order4;genus4;    62    genus    119


Dataset 1
*********
    - databank.fasta:
>A2 domain1;sub-domain1;class1;order1;sub-order1;family1;genus1;sp1
A
>B2 domain1;sub-domain1;class2;sub-class1;unidentified;family2;sp2
A
>C2 domain1;phylum1;class3;sub-order2;family3;genus2;sub-genus1;sp3
A
>D2 domain1;phylum1;sub-class2;order3;family4;genus3;sp4
A
>E2 domain1;phylum1;sub-class3;sub-family1;sp5
A
>F2 domain2;class4;sub-class4;order4;sub-order3;genus4;sp6
A

    - databank.tax:
domain1;    1    domain    119
domain1;sub-domain1;    2    sub-domain    119
domain1;sub-domain1;class1;    3    class    119
domain1;sub-domain1;class1;order1;    4    order    119
domain1;sub-domain1;class1;order1;sub-order1;    5    sub-order    119
domain1;sub-domain1;class1;order1;sub-order1;family1;    6    family    119
domain1;sub-domain1;class1;order1;sub-order1;family1;genus1;    7    genus    119
domain1;sub-domain1;class2;    8    class    119
domain1;sub-domain1;class2;sub-class1;    9    sub-class    119
domain1;sub-domain1;class2;sub-class1;unidentified;    10    order    119
domain1;sub-domain1;class2;sub-class1;unidentified;family2;    11    family    119
domain1;phylum1;    12    phylum    119
domain1;phylum1;class3;    13    class    119
domain1;phylum1;class3;sub-order2;    14    sub-order 119
domain1;phylum1;class3;sub-order2;family3;    15    family    119
domain1;phylum1;class3;sub-order2;family3;genus2;    16    genus    119
domain1;phylum1;class3;sub-order2;family3;genus2;sub-genus1;    17    sub-genus    119
domain1;phylum1;sub-class2;    18    sub-class    119
domain1;phylum1;sub-class2;order3;    19    order    119
domain1;phylum1;sub-class2;order3;family4;    20    family    119
domain1;phylum1;sub-class2;order3;family4;genus3;    21    genus    119
domain1;phylum1;sub-class3;    22    sub-class    119
domain1;phylum1;sub-class3;sub-family1;    23    sub-family    119
domain2;    24    domain    119
domain2;class4;    25    class    119
domain2;class4;sub-class4;    26    sub-class    119
domain2;class4;sub-class4;order4;    27    order    119
domain2;class4;sub-class4;order4;sub-order3;    28    sub-order    119
domain2;class4;sub-class4;order4;sub-order3;genus4;    29    genus    119
"""