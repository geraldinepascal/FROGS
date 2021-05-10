#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import os
import sys
import argparse
import gzip
import json
import re
import ete3 as ete
import inspect

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

PRO_DIR = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/prokaryotic/pro_ref/')
ITS_DIR = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_ITS/')


from frogsUtils import *
from frogsBiom import BiomIO
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES 
#
##################################################################################################################################################

class PlaceSeqs(Cmd):
	"""
	@summary: place studies sequences (i.e. OTUs) into a reference tree
	"""
	def __init__(self, in_fasta, out_tree, placement_tool, categorie):
		"""
		@param in_fasta: [str] Path to input fasta file.
		@param out_tree: [str] Path to output resulting tree file.
		@param placement_tool: [str] Placement tool to use (epa-ng or sepp).
		@param ref_dir: [str] Directory containing reference sequence files.
		"""
		if categorie == "16S":
			categorie = PRO_DIR
		elif categorie == "ITS":
			categorie = ITS_DIR

		Cmd.__init__(self,
		'place_seqs.py',
		'place OTU on reference tree.',
		'--study_fasta ' + in_fasta + ' --out_tree ' + out_tree + ' --placement_tool ' + placement_tool + " --ref_dir " + categorie,
		'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

class FindClosestsRefSequences(Cmd):
	'''
	@summary: find OTUs closests reference sequences into a reference tree.
	'''
	def __init__(self, in_tree, in_biom, out_summary):
		'''
		@param in_tree: [str] Path to place_seqs.py output tree.
		@param in_biom: [str] Path to BIOM input file.
		@summary out_summary: [str] Path to output summary file.
		'''
		Cmd.__init__(self,
			'find_closest_ref_sequence.py',
			'find OTUs closests reference sequences into a reference tree.',
			'--tree_file ' + in_tree + ' --biom_file ' + in_biom + ' --output ' + out_summary,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()	

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def convert_fasta(in_fasta, out_fasta):
	"""
	@summary: Change fasta headers to be compatible with picrust2
	"""
	FH_input = FastaIO(in_fasta)
	FH_output = FastaIO(out_fasta,"wt" )
	for record in FH_input:
		record.id = record.id
		record.description = None
		FH_output.write(record)
	FH_output.close()

def excluded_sequence(tree_file, in_fasta, excluded):
	"""
	@summary: Returns the excluded sequence.
	@param fasta_file: [str] Path to the fasta file to process.
	@param tree_file: [str] Path to the tree file to process.
	@return: [int] The file of no aligned sequence.
	"""
	file = open(tree_file, "r")
	line = file.readline()
	list_cluster = re.findall("(Cluster_[0-9]+)", line)
	file.close()

	FH_input = FastaIO(in_fasta)
	excluded = open(excluded, "wt")

	for record in FH_input:
		if record.id not in list_cluster:
			excluded.write(rrecord.id)
	FH_input.close()
	excluded.close()


def remove_excluded_fasta( in_fasta, out_fasta, excluded_seqs):
	excluded = []
	excluded_file = open(excluded_seqs)
	if excluded_file is not None:
		for li in excluded_file:
			excluded.append(li.strip())

	FH_input = FastaIO( in_fasta )
	FH_output = FastaIO( out_fasta, "wt")
	for record in FH_input:
		real_id = record.id.split()[0]

		if real_id not in excluded:
			FH_output.write(record)
	FH_input.close()
	FH_output.close()

def remove_observations(input_biom, output_biom, excluded_seqs):
	"""
	@summary: Removes the specified list of observations.
	@param removed_observations: [list] The names of the observations to remove.
	@param input_biom: [str] The path to the input BIOM.
	@param output_biom: [str] The path to the output BIOM.
	"""
	excluded = []
	excluded_file = open(excluded_seqs)
	if excluded_file is not None:
		for li in excluded_file:
			excluded.append(li.strip())

	biom = BiomIO.from_json( input_biom )
	biom.remove_observations( excluded )
	BiomIO.write( output_biom, biom )


def write_summary(summary_file, in_fasta, align_out, biomfile, closest_ref_files):
	"""
	@summary: Writes the process summary in one html file.
	@param summary_file: [str] path to the output html file.
	@param align_out: [str] path to the fasta file of unaligned OTU
	@param biomfile: [str] path to the input BIOM file.
	@param treefile: [str] path to the Newick file.
	"""
	# to summary OTUs number && abundances number			   
	summary_info = {
	   'nb_kept' : 0,
	   'nb_removed' : 0,
	   'abundance_kept' : 0,
	   'abundance_removed' : 0	   
	}
	number_otu_all = 0
	number_abundance_all = 0
	# to detail removed OTU
	details_categories =["Blast_taxonomy","Closest_ref_ID","Closest_ref_name","Closest_ref_taxonomy","Closest_ref_distance"]
	infos_otus = list()
	biom=BiomIO.from_json(biomfile)
	list_otu_all = []
	# record nb OTU and abundance
	for otu in FastaIO(in_fasta):
		list_otu_all.append(otu.id)
		number_otu_all +=1
		number_abundance_all += biom.get_observation_count(otu.id)

	closest_ref = open(closest_ref_files)
	for li in closest_ref:
		li = li.strip().split('\t')
		if "Closest_ref_ID" not in li:
			infos_otus.append({
				'name': li[0],
				'data': list(map(str,li[1:]))
				})

	# record details about removed OTU
	if align_out is not None:
		for otu in FastaIO(align_out):
			summary_info['nb_removed'] +=1
			summary_info['abundance_removed'] += biom.get_observation_count(otu.id)
	
	summary_info['nb_kept'] = number_otu_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']


	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep1_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )


	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categories) )
		elif "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_otus) )
		elif "###REMOVE_DATA###" in line:
			line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
		FH_summary_out.write( line )

	FH_summary_out.close()
	FH_summary_tpl.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="place studies sequences (i.e. OTUs) into a reference tree.")
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-i', '--input_fasta', required=True, help="Input fasta file of unaligned studies sequences")
	group_input.add_argument('-b', '--input_biom', required=True, help='Biom file.')
	group_input.add_argument('-c', '--categorie',choices=['16S', 'ITS'], default='16S', help='Specifies which categorie 16S or ITS')
	group_input.add_argument('-p', '--placement_tool', default='epa-ng', help='Placement tool to use when placing sequences into reference tree. One of "epa-ng" or "sepp" must be input')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--out_tree', default='out.tree', help='Tree output with insert sequences (format: newick).')
	group_output.add_argument('-e', '--excluded', default='excluded.fasta', help='List of sequences not inserted in the tree.')
	group_output.add_argument('-s', '--insert_sequences', default='place_seqs.fasta', help='sequences file without non insert sequences. (format: FASTA). [Default: %(default)s]')
	group_output.add_argument('-m', '--insert_biom', default='place_seqs.biom', help='abundance file without chimera (format: BIOM)')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.out_tree)[0])

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		Logger.static_write(args.log_file,'\n# Cleaning fasta headers\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
		tmp_fasta = tmp_files.add('cleaned.fasta')
		convert_fasta(args.input_fasta,tmp_fasta)

		if args.categorie == "ITS":

			gz_ref_fasta = os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_ITS/fungi_ITS.fna.gz')

			if os.path.exists(gz_ref_fasta):

				input_ref = gzip.GzipFile(gz_ref_fasta, 'rb')
				f = input_ref.read()
				input_ref.close()
				output = open(os.path.join(os.path.dirname(os.__file__),'site-packages/picrust2/default_files/fungi/fungi_ITS/fungi_ITS.fna'), 'wb')
				output.write(f)
				output.close()
				os.remove(gz_ref_fasta)


		try:
			PlaceSeqs(tmp_fasta, args.out_tree, args.placement_tool, args.categorie).submit(args.log_file)

		except subprocess.CalledProcessError:
			print('\n\n#ERROR : epa-ng running out of memory. Please use placement tool sepp instead ( -t sepp )')

		excluded_sequence(args.out_tree,args.input_fasta,args.excluded)
		remove_excluded_fasta(args.input_fasta, args.insert_sequences, args.excluded)
		remove_observations(args.input_biom, args.insert_biom, args.excluded)

		closest_ref_files = tmp_files.add( "closest_ref.tsv" )
		FindClosestsRefSequences(args.out_tree, args.input_biom, closest_ref_files).submit(args.log_file)
		write_summary(args.html, tmp_fasta, args.excluded, args.input_biom,closest_ref_files)

	finally:
		if not args.debug:
			tmp_files.deleteAll()
