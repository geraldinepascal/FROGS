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
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'dev'


import sys
import random
import argparse
import numpy as np
from frogsNode import Node
from frogsSequenceIO import FastaIO


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_subset( in_path, out_path, selected ):
    FH_in = FastaIO(in_path)
    FH_out = FastaIO(out_path, "w")
    for record in FH_in:
        if record.id in selected:
            FH_out.write(record)
    FH_in.close()
    FH_out.close()

def print_summary( tree ):
    max_depth = databank_tree.get_leaves()[0].get_depth()
    print "##########################################################\n" + \
          "# Representation\n" + \
          "#\n"
    print "Rank_depth\tNb_taxa\tNb_selected"
    for idx in range(max_depth):
        depth = idx + 1
        rank_nodes = databank_tree.get_descendants(depth)
        nb_node_selected = 0
        for node in rank_nodes:
            for leaf in node.get_leaves():
                if leaf.metadata["selected"]:
                    nb_node_selected += 1
                    break
        print depth, len(rank_nodes), nb_node_selected
    print ""

    print "##########################################################\n" + \
          "# Distribution\n" + \
          "#\n"
    distrib_data = dict()
    for idx in range(max_depth -1):
        depth = idx + 1
        nb_selected = list()
        rank_nodes = databank_tree.get_descendants(depth)
        for node in rank_nodes:
            nb_children_selected = 0
            for child in node.get_children():
                for leaf in child.get_leaves():
                    if leaf.metadata["selected"]:
                        nb_children_selected += 1
                        break
            nb_selected.append(nb_children_selected)
        distrib_data[str(depth)] = {
            "all": distrib(nb_selected),
            "non-zero": distrib(filter(lambda a: a != 0, nb_selected))
        }
    print "Distribution in all nodes:"
    print "\t" + "\t".join(sorted(distrib_data))
    for field in ["min", "10/100", "25/100", "50/100", "50/100", "75/100", "90/100", "max"]:
        print field + ":\t" + "\t".join([str(distrib_data[depth]["all"][field]) for depth in sorted(distrib_data)])
    print "\nDistribution in represented nodes:"
    print "\t" + "\t".join(sorted(distrib_data))
    for field in ["min", "10/100", "25/100", "50/100", "50/100", "75/100", "90/100", "max"]:
        print field + ":\t" + "\t".join([str(distrib_data[depth]["non-zero"][field]) for depth in sorted(distrib_data)])

def distrib( data ):
    return {
        "min": min(data),
        "10/100": np.percentile(data, 10),
        "25/100": np.percentile(data, 25),
        "50/100": np.percentile(data, 50),
        "75/100": np.percentile(data, 75),
        "90/100": np.percentile(data, 90),
        "max": max(data)
    }

def ascending_walk(node, max_selection):
    selected_leaf = list()
    if max_selection > 0:
        if random.randint(1, args.climb_prob) != 1:
            log[-1].append("ascending")
            parent = node.get_parent()
            if parent is not None: # Node is not root
                brothers = parent.get_children()
                if len(brothers) > 1: # Node has brother(s)
                    if random.randint(1, args.neighbor_prob) != 1: # Brother recruitment
                        neighbors_leaves = list()
                        for brother in brothers:
                            if brother is not node:
                                for leaf in brother.get_leaves():
                                    if not leaf.metadata["selected"]:
                                        neighbors_leaves.append( leaf )
                        if len(neighbors_leaves) > 0:
                            selected_idx = random.randint(1, len(neighbors_leaves)) -1
                            log[-1].append("neighbor_selection: " + neighbors_leaves[selected_idx].metadata["retained_seq_id"])
                            selected_leaf.append( neighbors_leaves[selected_idx] )
                            max_selection -= 1
                    # Go to parent
                    selected_leaf.extend( ascending_walk(parent, max_selection) )
    return selected_leaf


def descending_walk( node ):
    selected_leaf = None
    if not node.has_child():
        selected_leaf = node
    else:
        accessible_children = list()
        for child in node.get_children():
            for leaf in child.get_leaves():
                if not leaf.metadata["selected"]:
                    accessible_children.append(child)
                    break
        selected_idx = random.randint(1, len(accessible_children)) -1
        selected_leaf = descending_walk( accessible_children[selected_idx] )
    return selected_leaf


def rank_sampling( tree, rank ):
    selected_leaf = None
    accessible_leaves = list()
    for rank_node in tree.get_descendants(rank):
        rank_is_accessible = rank_node.has_child
        for leaf in rank_node.get_leaves():
            if leaf.metadata["selected"]:
                rank_is_accessible = False
                break
        if rank_is_accessible:
            accessible_leaves.extend( rank_node.get_leaves() )
    selected_idx = random.randint(1, len(accessible_leaves)) -1
    selected_leaf = accessible_leaves[selected_idx]
    return selected_leaf


def get_tree_from_fasta( in_fasta ):
    """
    @warning: The root node must be present
    """
    databank_tree = None
    FH_databank = FastaIO(in_fasta)
    for record in FH_databank:
        if record.description.endswith(";"):
            record.description = record.description[:-1]
        taxonomy = record.description.split(";")
        if databank_tree is None:
            databank_tree = Node(taxonomy[0])
        parent = databank_tree
        for rank_depth, taxa in enumerate(taxonomy[1:]):
            if not parent.has_child( taxa ):
                taxa_node = Node(taxa, parent)
                if (rank_depth+1) == (len(taxonomy)-1): # Current node is leaf
                    taxa_node.metadata["seq_ids"] = [record.id]
            else:
                if (rank_depth+1) == (len(taxonomy)-1): # Current node is leaf
                    taxa_node = parent.get_child(taxa)
                    taxa_node.metadata["seq_ids"].append(record.id)
            parent = parent.get_child(taxa)
    FH_databank.close()
    return databank_tree


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=('''Produce a random subset of species (with near and distant species) from a databank.''')
    )
    parser.add_argument( '-c', '--climb-prob', default=4, type=int, help='Porbability when you have selected a node that you go to the parent node (to eventualy select a neighbor node) [DEFAULT: 4]. Example: -c 4 you have 1/4 chance to stop ascension ; -c 5 you have 1/4 chance to stop ascension.' )
    parser.add_argument( '-n', '--neighbor-prob', default=4, type=int, help='Porbability when you have go to a parent to select a neighbor nodes in children [DEFAULT: 4]. Example: -c 4 you have 1/4 chance to skip selection at this level ; -c 5 you have 1/4 chance to skip selection at this level.' )
    parser.add_argument( '-e', '--expected-nb', default=100, type=int, help='Number of selected sequences.' )
    parser.add_argument( '-f', '--select-from', default="bottom", choices=['top', 'bottom', 'mix'], help='Select in "top", "bottom" or "mix". With top the less divided branch is favored ; with bottom the most divided branch is favored ; with mix the less divided and most divided are favored. [Default: bottom]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-d', '--databank', required=True, help='The reference databank (format : FASTA). Each sequence must have the same number of tacxonomy level and the header must have this format: "ID<TAB>TAX_LVL1;TAX_LVL2".' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output', required=True, help='The selected sequences (format : FASTA).' )
    args = parser.parse_args()

    log = list()

    # Build tree
    databank_tree = get_tree_from_fasta(args.databank)

    # Select one sequence ID by leaf
    for leaf in databank_tree.get_leaves():
        nb_seq_ids = len(leaf.metadata["seq_ids"])
        leaf.metadata["selected"] = False
        leaf.metadata["retained_seq_id"] = leaf.metadata["seq_ids"][0]
        if nb_seq_ids > 1:
            leaf.metadata["retained_seq_id"] = leaf.metadata["seq_ids"][random.randint(0, nb_seq_ids-1)]

    # Select leaves
    current_nb = 0
    next_from_top = (args.select_from == "top")
    nb_asc = 0
    selected_leaves_id = list()
    while args.expected_nb > current_nb:
        # Random selection
        current_selected_leaf = None
        if next_from_top:
            current_selected_leaf = descending_walk(databank_tree)
            log.append(["from_top_selection: " + current_selected_leaf.metadata["retained_seq_id"]])
        else:
            current_selected_leaf = rank_sampling( databank_tree, 6 )#################################################### Param
            log.append(["from_bottom_selection: " + current_selected_leaf.metadata["retained_seq_id"]])
            nb_asc += 1
        if args.select_from == "mix":
            if nb_asc == 2:
                nb_asc = 0
                next_from_top = True
            else:
                next_from_top = False
        current_selected_leaf.metadata["selected"] = True
        selected_leaves_id.append( current_selected_leaf.metadata["retained_seq_id"] )
        current_nb += 1

        # Neighbor selection
        current_selected_leaves = ascending_walk(current_selected_leaf, (args.expected_nb-current_nb))
        for leaf in current_selected_leaves:
            leaf.metadata["selected"] = True
            selected_leaves_id.append( leaf.metadata["retained_seq_id"] )
            current_nb += 1

    # Write selection
    write_subset(args.databank, args.output, selected_leaves_id)

    # Log
    for action in log:
        print action

    # Summary
    print_summary(databank_tree)
