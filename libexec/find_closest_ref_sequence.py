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

from frogsUtils import *
from frogsBiom import BiomIO
from frogsSequenceIO import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def find_clusters(tree):
	"""
	@summary: Returns the list of clusters insert into reference tree.
	@param tree: [str] Path to tree output from place_seqs.py.
	"""
	fi = open(tree,'r').readline()
	list_cluster = re.findall("(Cluster_[0-9]+)", fi)
	return list_cluster

def find_closest_ref_sequences(tree, biom_file, fasta_file, clusters, ref_file, output):
	"""
	@summary: Find closest reference sequence in the tree for every cluster.
	@param tree: [str] Path to tree output from place_seqs.py.
	@param biom_file: [str] path to BIOM input file.
	@param fasta_file: [str] path to fasta input file.
	@param clusters: [list] clusters insert in tree (find_clusters output).
	@ref_file: [str] path to reference map file in order to have taxonomies informations.
    """
	biom=BiomIO.from_json(biom_file)
	if args.category == '16S':
		ref_file = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "tools/FPStep1/data/JGI_ID_to_taxonomy.txt"))
		picrust_aln = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/prokaryotic/pro_ref/pro_ref.fna')

	elif args.category == 'ITS':
		ref_file = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "tools/FPStep1/data/JGI_ID_ITS_to_taxonomy.txt"))
		picrust_aln = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_ITS/fungi_ITS.fna')
	
	ref = open(ref_file,'r').readlines()
	ID_to_taxo = {}
	for li in ref[1:]:
		li = li.strip().split('\t')
		#ID of reference sequence to [sequence name, taxonomy]
		ID_to_taxo[li[0]] = [li[1],li[2]]

	ref_seqs = {}
	FH_picrust2_aln = FastaIO( picrust_aln )
	for record in FH_picrust2_aln:
		record.id = record.id.split('-')[0]
		record.string = record.string.replace('-','').upper()
		ref_seqs[record.id] = record.string

	cluster_to_seq = {}
	FH_input = FastaIO( fasta_file )
	for record in FH_input:
		cluster_to_seq[record.id] = record.string

	FH_out = open(output,'wt')
	FH_out.write('Cluster\tTaxonomy\tClosest_ref_ID\tClosest_ref_name\tClosest_ref_taxonomy\tClosest_ref_distance\tComment')

	t=ete.Tree(tree)

	for cluster in clusters:
		FH_out.write(cluster+'\t'+";".join(biom.get_observation_metadata(cluster)["blast_taxonomy"])+'\t')

		node = t.search_nodes(name=cluster)[0]
		#find distances from cluster to every reference sequences is sister group.
		for sister_group in node.get_sisters():
			leaf_to_dist = {}
			for leaf in sister_group.get_leaves():
				leaf.name = leaf.name.replace('-cluster','')
				# if sequence in sister group is not another cluster
				if leaf.name not in clusters:
					leaf_to_dist[leaf.name] = t.get_distance(leaf,cluster)

			best_leaf = sorted(leaf_to_dist, key=leaf_to_dist.get)[0]
			if best_leaf in ID_to_taxo:
				#cleaning leaf name
				best_leaf = best_leaf.split('-')[0]
				comment = "/"
				genus_picrust = ID_to_taxo[best_leaf][1].split(';')[-2]
				species_picrust = " ".join(ID_to_taxo[best_leaf][1].split(';')[-1].split(' ')[0:2])
				genus_frogs = biom.get_observation_metadata(cluster)["blast_taxonomy"][-2]
				species_frogs = biom.get_observation_metadata(cluster)["blast_taxonomy"][-1]
				if genus_picrust == genus_frogs and species_picrust == species_frogs:
					comment = "identical taxonomy"
				# if 100% identity on OTU length against reference
				if cluster_to_seq[cluster] in ref_seqs[best_leaf]:
					if comment == "/":
						comment = "identical sequence"
					else:
						comment+=";identical sequence"
				FH_out.write(best_leaf+'\t'+ID_to_taxo[best_leaf][0]+'\t'+ID_to_taxo[best_leaf][1]+'\t'+str(leaf_to_dist[best_leaf])+'\t'+str(comment)+'\n')

			else:
				FH_out.write(' \t \t \t \t \n')
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
	group_input.add_argument('-f', '--fasta_file', required=True, help='Input fasta file.')
	group_input.add_argument('-b', '--biom_file', required=True, help='Input biom file.')
	group_input.add_argument('-c', '--category', choices=['16S', 'ITS'], default='16S', help='Specifies which category 16S or ITS')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output = parser.add_argument('-o', '--output', default='closests_ref_sequences.txt')
	group_output = parser.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

	args = parser.parse_args()
	prevent_shell_injections(args)

	clusters = find_clusters(args.tree_file)

	find_closest_ref_sequences(args.tree_file, args.biom_file, args.fasta_file, clusters, args.category, args.output)
	