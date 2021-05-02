#!/usr/bin/env python3

__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os, sys
import argparse
import re
import ete3 as ete

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
os.environ['PATH'] = CURRENT_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

REF_FILE = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "tools/FPStep1/data/JGI_ID_to_taxonomy.txt"))

from frogsUtils import *
from frogsBiom import BiomIO
from frogsSequenceIO import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def find_clusters(tree):

	fi = open(tree,'r').readline()

	list_cluster = re.findall("(Cluster_[0-9]+)", fi)

	return list_cluster

def find_closests_ref_sequences(tree, biom_file, clusters, ref_file, output):

	biom=BiomIO.from_json(biom_file)
	ref = open(ref_file,'r').readlines()
	ID_to_taxo = {}

	for li in ref[1:]:
		li = li.strip().split('\t')
		ID_to_taxo[li[0]] = li[1]

	FH_out = open(output,'wt')
	FH_out.write('Cluster\tBlast_taxonomy\tclosest_ref_ID\tclosest_ref_taxonomy\tclosest_ref_distance\n')

	t=ete.Tree(tree)

	for cluster in clusters:
		FH_out.write(cluster+'\t'+";".join(biom.get_observation_metadata(cluster)["blast_taxonomy"])+'\t')

		node = t.search_nodes(name=cluster)[0]

		for sister_group in node.get_sisters():
			leaf_to_dist = {}
			for leaf in sister_group.get_leaves():
				leaf_to_dist[leaf.name] = sister_group.get_distance(leaf)

			for leaf in sorted(leaf_to_dist, key=leaf_to_dist.get):
				if leaf not in clusters and leaf in ID_to_taxo:
					#cleaning leaf name
					leaf = leaf.split('-')[0]

					FH_out.write(leaf+'\t'+ID_to_taxo[leaf]+'\t'+str(leaf_to_dist[leaf])+'\n')
					break

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="Find OTUs closests references sequences in tree")
	parser.add_argument( '-v', '--version', action='version', version=__version__ )

	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-t', '--tree_file', required=True, help='Tree file (output of place_seqs.py')
	group_input.add_argument('-b', '--biom_file', required=True, help='Biom file.')
	group_input.add_argument('-r', '--ref_sequences', default=REF_FILE, help='JGI reference file (JGI ID to JGI taxonomy)')

	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output = parser.add_argument('-o', '--output', default='closests_ref_sequences.txt')
	group_output = parser.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

	args = parser.parse_args()
	prevent_shell_injections(args)

	clusters = find_clusters(args.tree_file)

	find_closests_ref_sequences(args.tree_file, args.biom_file, clusters, args.ref_sequences, args.output)
	