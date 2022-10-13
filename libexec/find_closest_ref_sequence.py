#!/usr/bin/env python3
#
# Copyright (C) 2022 INRA
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
__author__ = 'Vincent Darbot INRAE - GENPHYSE'
__copyright__ = 'Copyright (C) 2022 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'dev'

from Bio import Align
import ete3 as ete
import argparse
import os, sys
import random
import gzip
import re

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
	if taxo_frogs == "NA":
		return False
	for i in range(len(taxo_frogs.split(';'))):
		if taxo_frogs.split(';')[i].lower() != taxo_picrust.split(';')[i].lower():
			return False
	return True

def find_lowest_same_taxo_rank(taxo_frogs, taxo_picrust, hierarchy = ["Up to Kingdom","Up to Phylum","Up to Class","Up to Order","Up to Family","Up to Genus","Up to Species"]):
	'''
	@summary: find lowest identical taxonomic rank between frogs and picrust2 taxonomies
	return rank_level and rank_name
	'''
	if taxo_frogs == "NA":
		return "/"
	taxo_frogs = [taxo_frogs.split(';')[i].lower() for i in range(len(taxo_frogs.split(';')))]
	taxo_picrust = [taxo_picrust.split(';')[i].lower() for i in range(len(taxo_picrust.split(';')))]
	for i in range(len(hierarchy)-1, -1, -1):
		if taxo_frogs[i] == taxo_picrust[i]:
			return i, hierarchy[i]
	return 7, "/"

def iter_sample_fast(iterator, samplesize):
    results = []
    
    # Fill in the first samplesize elements:
    for _ in range(samplesize):
        results.append(iterator.__next__())
    random.shuffle(results)  # Randomize their positions
    
    return results[:samplesize]

def run_megablast(query, subject):
	'''
	@summary: Run megablast to calculate %identity between frogs sequence \
	and picrust2 closest ref.
	return blast results: number of alignments, % of identity, % of aligment coverage, blast score
	'''
	## default paramaters
	aligner = Align.PairwiseAligner()
	aligner.match_score = 1.0
	aligner.mismatch_score = -2.0
	aligner.gap_score = -2.5
	aligner.mode = 'local'
	######
	alignments = aligner.align(query, subject)
	n = str(len(alignments)) if len(alignments) <=500 else ">500"
	scores = list()
	covs = list()
	identities = list()
	# check up to the first 500 first alignment
	if len(alignments) > 500:
		alignments =  iter_sample_fast(alignments, 500)
	for aln in alignments:
		if not aln.score in scores:
			scores.append(aln.score)
		coded_aln = aln.format().split('\n')[1].strip()
		n_match = coded_aln.count('|')
		n_gap = coded_aln.count('-')
		n_mismatch = coded_aln.count('.')

		qstart = aln.path[0][0]
		qend = aln.path[-1][0]

		cov = round(abs(qend-qstart)*100.0/len(query),2)
		if not cov in covs:
			covs.append(cov)
		identity = round(n_match * 100.0 / abs(qend-qstart),2)
		if not identity in identities:
			identities.append(identity)

	id = str(min(identities)) + " - " + str(max(identities)) if len(identities) > 1 else str(identities[0])
	cov = str(min(covs)) + " - " + str(max(covs)) if len(covs) > 1 else str(covs[0])
	score = str(min(scores)) + " - " + str(max(scores)) if len(scores) > 1 else str(scores[0])

	res = {
		"n_aln" : n,
		"id" : id,
		"cov" : cov,
		"score" : score
	}

	return res

def check_ref_files(tree_file, biom_file, biom_path, fasta_file, ref_aln, output ):
	'''
	@param tree: [str] Path to tree output from place_seqs.py.
	@param biom_file: [str] path to BIOM input file.
	@param fasta_file: [str] path to fasta input file.
	@ref_file: [str] path to reference map file in order to have taxonomies informations.
	'''
	biom=BiomIO.from_json(biom_file)
	ref_file = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "frogsfunc_suppdata/JGI_ID_to_taxonomy.txt.gz"))
	picrust_aln = ref_aln
	
	FH_in = gzip.open(ref_file,'rt')
	ID_to_taxo = dict()
	for line in FH_in:
		line = line.strip().split('\t')
		#ID of reference sequence to [sequence name, taxonomy]
		ID_to_taxo[line[0]] = [line[1],line[2]]
	FH_in.close()

	ref_seqs = dict()
	#Assign each picrust ref sequence name to its sequence
	FH_picrust2_aln = FastaIO( picrust_aln )
	for record in FH_picrust2_aln:
		record.id = record.id.split('-')[0]
		record.string = record.string.replace('-','').upper()
		ref_seqs[record.id] = record.string
	FH_picrust2_aln.close()

	cluster_to_seq = dict()
	#Assign each FROGS cluster sequence name to its sequence
	FH_input = FastaIO( fasta_file )
	for record in FH_input:
		cluster_to_seq[record.id] = record.string
	FH_input.close()
	# ete3 input tree file
	tree=ete.Tree(tree_file)

	clusters = list()
	for leaf in tree:
		if leaf.name.split('-cluster')[0] not in ID_to_taxo.keys():
			clusters.append(leaf.name)
	return [tree, biom, biom_path, ID_to_taxo, ref_seqs, cluster_to_seq, output, clusters]

def find_closest_ref_sequences(tree, biom, biom_path, ID_to_taxo, ref_seqs, cluster_to_seq, output, clusters):
	"""
	@summary: find each closest picrust ref sequence from FROGS cluster, from FPStep1 tree output file.
	@param tree: [str]  Tree as input for ete3.
	@param biom: [str] Biom read from jason file.
	@param ID_to_taxo: [dico] Associate each picrust2 ref sequence to its taxonomy.
	@param ref_seqs: [dico] Assign each picrust ref sequence name to its sequence.
	@param clusters_to_seq: [list] clusters insert in tree and sequences associated(find_clusters output).
	@param output: [str] path to tmp output file in order to write frogs and picrust2 taxonomic comparaisons.
	"""
	max_nsti = 0
	FH_out = open(output,'wt')
	header = "\t".join(["#Cluster","Nb sequences", "FROGS Taxonomy",\
		"PICRUSt2 closest ID","PICRUSt2 closest reference name","PICRUSt2 closest taxonomy",\
		"NSTI", "NSTI Confidence" ,"FROGS and PICRUSt2 lowest same taxonomic rank",\
		 "Comment", "Cluster sequence", "PICRUSt2 closest reference sequence",\
		"n_aln", "%id", "%cov", "score"])
	FH_out.write(header+"\n")

	for observation_name in biom.get_observations_names():
		count = str(biom.get_observation_count(observation_name))
		node = tree.search_nodes(name=observation_name)[0]
		frogs_taxo = 'NA'
		affis_picrust = 'NA'
		ref_leaf_id = ''
		ref_leaf_taxo = ''
		lowest_same_rank = '/'
		comment = '/'
		blast_n_aln = ''
		blast_id = ''
		blast_cov = ''
		blast_score = ''

		leaf_to_dist = dict()
		for sister_group in node.get_sisters():
			for leaf in sister_group.get_leaves():
				leaf.name = leaf.name.replace('-cluster','')
				# if sequence in sister group is not another cluster
				if leaf.name not in clusters:
					leaf_to_dist[leaf.name] = tree.get_distance(leaf, observation_name)
		best_leaf = sorted(leaf_to_dist, key=leaf_to_dist.get)[0]

		if best_leaf in ID_to_taxo:
			ref_leaf_id = ID_to_taxo[best_leaf][0]
			ref_leaf_taxo = ID_to_taxo[best_leaf][1]
			affis_picrust = ID_to_taxo[best_leaf][1].replace(' ','_')
			rank_level = 7
			found_same_taxo = False

			if ( not biom.get_observation_metadata(observation_name)["blast_taxonomy"] is None or not len(biom.get_observation_metadata(observation_name)["blast_taxonomy"]) == 0 ) and not found_same_taxo:
				for affi in biom.get_observation_metadata(observation_name)["blast_affiliations"]:
					cur_frogs_taxo = ";".join(affi["taxonomy"])
					if '__' in cur_frogs_taxo:
						cur_frogs_taxo = ";".join(["".join(af.split('__')[1:]) for af in cur_frogs_taxo.split(';')])
					clean_cur_frogs_taxo = cur_frogs_taxo.replace(' ','_')
					cur_rank_level, cur_lowest_same_rank = find_lowest_same_taxo_rank(clean_cur_frogs_taxo, affis_picrust)

					if cur_rank_level < rank_level:
						rank_level = cur_rank_level
						lowest_same_rank = cur_lowest_same_rank
						frogs_taxo = cur_frogs_taxo

						if is_same_taxonomies(frogs_taxo, affis_picrust):
							found_same_taxo = True
							comment = "identical taxonomy"

		biom.add_metadata(observation_name, "picrust2_affiliations", affis_picrust, "observation", erase_warning = False)
		
		if cluster_to_seq[observation_name] in ref_seqs[best_leaf]:
			if comment == "/":
				comment = "identical sequence"
			else:
				comment+=";identical sequence"

		blast = run_megablast(cluster_to_seq[observation_name], ref_seqs[best_leaf])
		blast_n_aln, blast_id, blast_cov, blast_score = blast['n_aln'], blast['id'], blast['cov'], blast['score']

		confidence = "To exclude"
		if rounding(leaf_to_dist[best_leaf]) >= 1 and rounding(leaf_to_dist[best_leaf]) < 2:
			confidence = "Bad"
		elif rounding(leaf_to_dist[best_leaf]) >= 0.5 and rounding(leaf_to_dist[best_leaf]) < 1:
			confidence = "Medium"
		elif rounding(leaf_to_dist[best_leaf]) < 0.5:
			confidence = "Good"

		if rounding(leaf_to_dist[best_leaf]) > max_nsti:
			max_nsti = rounding(leaf_to_dist[best_leaf])

		FH_out.write("\t".join([observation_name, count, frogs_taxo, best_leaf,\
		ref_leaf_id, ref_leaf_taxo, str(rounding(leaf_to_dist[best_leaf])),\
		confidence, lowest_same_rank, comment, cluster_to_seq[observation_name], ref_seqs[best_leaf],\
		blast_n_aln, blast_id, blast_cov, blast_score])+'\n')
	BiomIO.write(biom_path, biom)
	return max_nsti

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
	group_input.add_argument('-t', '--input-tree', required=True, help='Tree file (output of place_seqs.py')
	group_input.add_argument('-f', '--input-fasta', required=True, help='Input fasta of sequences included in PICRUSt2 reference tree.')
	group_input.add_argument('-b', '--input-biom', required=True, help='Input biom file of sequences included in PICRUSt2 reference tree.')
	group_input.add_argument('-r', '--ref-aln', required=True, help='Alignment of reference sequences used in FPStep1 in order to execute place_seqs.py (ie $PICRUST_PATH/frogsfunc_suppdata/fungi/fungi_ITS/')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output = parser.add_argument('-o', '--output-tsv', default='FPStep1_closests_ref_sequences.tsv', help='Informations about clusters and picrust2 closest reference from cluster sequences (identifiants, taxonomies, phylogenetic distance from reference, nucleotidics sequences')
	group_output = parser.add_argument('-e', '--output-biom', default='FPStep1.biom', help='Biom file without non insert sequences. (format: BIOM) [Default: %(default)s]')
	group_output = parser.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')

	args = parser.parse_args()
	Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
	prevent_shell_injections(args)

	inputs = check_ref_files(args.input_tree, args.input_biom, args.output_biom, args.input_fasta, args.ref_aln, args.output_tsv )
	tree_formatted = inputs[0]
	reference_sequences = inputs[3]

	clusters = inputs[-1]

	max_nsti = find_closest_ref_sequences(*inputs)
	Logger.static_write(args.log_file, "#Max NSTI: "+str(max_nsti))
	
