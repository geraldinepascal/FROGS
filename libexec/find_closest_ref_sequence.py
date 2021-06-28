#!/usr/bin/env python3

__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'dev'

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

def rounding(nb):
	'''
	@summary: Rounding numbers decimal 
	'''
	if re.search("^[0-9]{1}[.][0-9]+e",str(nb)):
		start = re.compile("[0-9][.][0-9]{1,2}")
		end = re.compile("e-[0-9]+")
		return float("".join(start.findall(str(nb))+end.findall(str(nb))))

	elif re.search("[0][.][0]+",str(nb)):
		motif = re.compile("[0][.][0]+[0-9]{2}")
		return float("".join(motif.findall(str(nb))))

	elif re.search("[0][.][1-9]+",str(nb)):
		return(round(nb,2))

	else:
		return nb

def is_same_taxonomies(taxo_frogs, taxo_picrust):
	'''
	@summary: compare if frogs and picrust taxonomies are egal are not
	@note: taxo inputs must be on this format: Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium_antarcticum
	'''
	for i in range(len(taxo_frogs.split(';'))):

		if taxo_frogs.split(';')[i].lower() != taxo_picrust.split(';')[i].lower():
			return False
	return True

def check_ref_files(tree_file, biom_file, multi_affi_file, fasta_file, category, output ):
	'''
	@param tree: [str] Path to tree output from place_seqs.py.
	@param biom_file: [str] path to BIOM input file.
	@param multi_affi: [str] path to multi-affiliations from biom input file. Run multiAffiFromBiom.py to generate this input.
	@param fasta_file: [str] path to fasta input file.
	@ref_file: [str] path to reference map file in order to have taxonomies informations.
	'''
	biom=BiomIO.from_json(biom_file)
	if category == '16S':
		ref_file = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "tools/FPStep1/data/JGI_ID_to_taxonomy.txt"))
		picrust_aln = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/prokaryotic/pro_ref/pro_ref.fna')

	elif category == 'ITS' or category == '18S':
		ref_file = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "tools/FPStep1/data/new_JGI_ITS_to_taxonomy.txt"))
		if category == 'ITS':
			picrust_aln = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_ITS/fungi_ITS.fna')
		else:
			picrust_aln = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_18S/fungi_ITS.fna')
	
	ref = open(ref_file,'r').readlines()
	ID_to_taxo = {}
	for li in ref[1:]:
		li = li.strip().split('\t')
		#ID of reference sequence to [sequence name, taxonomy]
		ID_to_taxo[li[0]] = [li[1],li[2]]

	ref_seqs = {}
	#Assign each picrust ref sequence name to its sequence
	FH_picrust2_aln = FastaIO( picrust_aln )
	for record in FH_picrust2_aln:
		record.id = record.id.split('-')[0]
		record.string = record.string.replace('-','').upper()
		ref_seqs[record.id] = record.string

	cluster_to_seq = {}
	#Assign each FROGS cluster sequence name to its sequence
	FH_input = FastaIO( fasta_file )
	for record in FH_input:
		cluster_to_seq[record.id] = record.string

	cluster_to_multiaffi = {}
	# Find each affiliation of FROGS multi-affiliations cluster 
	multi_affiliations = open(multi_affi_file,'r').readlines()
	for li in multi_affiliations:
		li = li.strip().split('\t')
		cluster = li[0]
		affi = li[1]
		if cluster not in cluster_to_multiaffi:
			cluster_to_multiaffi[cluster] = [affi]
		else:
			cluster_to_multiaffi[cluster].append(affi)

	# ete3 input tree file
	tree=ete.Tree(tree_file)

	return [tree, biom, cluster_to_multiaffi, ID_to_taxo, ref_seqs, cluster_to_seq, output]

def find_closest_ref_sequences(tree, biom, cluster_to_multiaffi, ID_to_taxo, ref_seqs, cluster_to_seq, output):
	"""
	@summary: find each closest picrust ref sequence from FROGS cluster, from FPStep1 tree output file.
	@param tree: Tree as input for ete3.
	@param biom: Biom read from jason file.
	@param biom cluster_to_multi_affi: [dico] Associate each FROGS cluster to its multi-affiliations.
	@param ID_to_taxo: [dico] Associate each picrust2 ref sequence to its taxonomy.
	@param ref_seqs: [dico] Assign each picrust ref sequence name to its sequence.
	@param clusters: [list] clusters insert in tree (find_clusters output).
	@ref_file: [str] path to reference map file in order to have taxonomies informations.
    """

	FH_out = open(output,'wt')
	FH_out.write('Cluster\tTaxonomy\tClosest_ref_ID\tClosest_ref_name\tClosest_ref_taxonomy\tClosest_ref_distance\tComment')

	for cluster in clusters:
		FH_out.write(cluster+'\t'+";".join(biom.get_observation_metadata(cluster)["blast_taxonomy"])+'\t')

		node = tree.search_nodes(name=cluster)[0]
		#find distances from cluster to every reference sequences is sister group.
		for sister_group in node.get_sisters():
			leaf_to_dist = {}
			for leaf in sister_group.get_leaves():
				leaf.name = leaf.name.replace('-cluster','')
				# if sequence in sister group is not another cluster
				if leaf.name not in clusters:
					leaf_to_dist[leaf.name] = tree.get_distance(leaf,cluster)

			best_leaf = sorted(leaf_to_dist, key=leaf_to_dist.get)[0]

			if best_leaf in ID_to_taxo:
				#cleaning leaf name
				best_leaf = best_leaf.split('-')[0]
				comment = "/"
				affis_picrust = ID_to_taxo[best_leaf][1].replace(' ','_')

				if cluster in cluster_to_multiaffi:
					for affi in cluster_to_multiaffi[cluster]:
						#formate FROGS taxonomy when when it's k__Fungi. k__Fungi --> Fungi
						if '__' in affi:
							affis_frogs = ";".join(["".join(af.split('__')[1:]) for af in affi.split(';')])
						else:
							affis_frogs = affi.replace(' ','_')

						if is_same_taxonomies(affis_frogs, affis_picrust):
							comment = "identical taxonomy"
							break
				else:
					affis_frogs = ";".join(biom.get_observation_metadata(cluster)["blast_taxonomy"])
					if '__' in affi:
						affis_frogs = ";".join(["_".join(af.split('__')[1:]) for af in affis_frogs.split(';')])
					else:
						affis_frogs = affis_frogs.replace(' ','_')

					if is_same_taxonomies(affis_frogs, affis_picrust):
						comment = "identical taxonomy"
						
				if cluster_to_seq[cluster] in ref_seqs[best_leaf]:
					if comment == "/":
						comment = "identical sequence"
					else:
						comment+=";identical sequence"
				FH_out.write(best_leaf+'\t'+ID_to_taxo[best_leaf][0]+'\t'+ID_to_taxo[best_leaf][1]+'\t'+str(rounding(leaf_to_dist[best_leaf]))+'\t'+str(comment)+'\n')
				
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
	group_input.add_argument('-m', '--multi_affi', required=True, help='Multi-affiliations from biom input file. Run multiAffiFromBiom.py to generate this input.')
	group_input.add_argument('-c', '--category', choices=['16S', 'ITS', '18S'], default='16S', help='Specifies which category 16S, ITS, 18S')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output = parser.add_argument('-o', '--output', default='closests_ref_sequences.txt')
	group_output = parser.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

	args = parser.parse_args()
	prevent_shell_injections(args)

	clusters = find_clusters(args.tree_file)

	inputs = check_ref_files(args.tree_file, args.biom_file, args.multi_affi, args.fasta_file, args.category, args.output )

	find_closest_ref_sequences(*inputs)
	