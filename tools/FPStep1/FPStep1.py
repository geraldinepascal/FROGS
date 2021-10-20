#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard & Vincent Darbot & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import os
import re
import sys
import gzip
import json
import inspect
import argparse
import ete3 as ete

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPAT
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

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
	def __init__(self, in_fasta, out_tree, placement_tool, ref_dir, min_align, log):
		"""
		@param in_fasta: [str] Path to input fasta file of unaligned cluster sequences.
		@param out_tree: [str] Path to output resulting tree file with insert clusters sequences.
		@param placement_tool: [str] Placement tool to use (epa-ng or sepp).
		@param ref_dir: [str] Directory containing reference sequence files.
		"""
		if ref_dir is None:
			opt = ''
		else:
			opt = ' --ref_dir ' + ref_dir

		Cmd.__init__(self,
		'place_seqs.py',
		'Place studies sequences (i.e. OTUs) on reference tree.',
		'--study_fasta ' + in_fasta + ' --out_tree ' + out_tree + ' --placement_tool ' + placement_tool + " --min_align " + str(min_align) + opt + " --verbose 2>> " + log,
		'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

class multiAffiFromBiom(Cmd):
	'''
    @summary: Extracts multi-affiliations data from a FROGS BIOM file, in order to find dentical sequences for FindClosestRefSequences step.
    @param in_biom: [Biom] The BIOM object.
    @param output_tsv: [str] Path to the output file (format : TSV).
    '''
	def __init__(self, in_biom, out_tsv, log):

		Cmd.__init__(self,
			'multiAffiFromBiom.py',
			'Extracts multi-affiliations data from a FROGS BIOM file.',
			'--input-file ' + in_biom + ' --output-file ' + out_tsv + " 2> " + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class FindClosestsRefSequences(Cmd):
	'''
	@summary: find OTUs closest reference sequences into a reference tree.
	'''
	def __init__(self, in_tree, in_biom, multi_affiliations, in_fasta, ref_aln, out_biom, out_summary):
		'''
		@param in_tree: [str] Path to resulting tree file with insert clusters sequences.(place_seqs.py output).
		@param in_biom: [str] Path to BIOM input file.
		@param category: [str] ITS or 16S
		@summary out_summary: [str] Path to output summary file.
		@note : Header of summary file:
		Cluster	Taxonomy	Closest_ref_ID	Closest_ref_name	Closest_ref_taxonomy	Closest_ref_distance
		'''
		Cmd.__init__(self,
			'find_closest_ref_sequence.py',
			'find OTUs closests reference sequences into a reference tree.',
			'--tree_file ' + in_tree + ' --biom_file ' + in_biom + ' --multi_affi ' + multi_affiliations + ' --fasta_file ' + in_fasta + ' --ref_aln ' + ref_aln + ' --out_biom ' + out_biom + ' --output ' + out_summary ,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class RemoveSeqsBiomFasta(Cmd):
	'''
	@summary: Create new biom and fasta file without not insert sequences in tree.
	'''
	def __init__(self, in_fasta, in_biom, out_fasta, out_biom, excluded_file):
		'''
		@param in_fasta: [str] Path to fasta input file.
		@param in_biom: [str] Path to BIOM input file.
		@param out_fasta: [str] Path to fasta output file.
		@param out_biom: [str] Path to BIOM output file.
		@param excluded_file: [str] Path to not insert sequences file (one Cluster ID per line).
		'''
		Cmd.__init__(self,
			'remove_seqs_biom_fasta.py',
			'remove not insert sequences in tree from fasta and biom file.',
			'--input_biom ' + in_biom + ' --input_fasta ' + in_fasta + ' --excluded_sequences ' + excluded_file + ' --output_biom ' + out_biom + " --output_fasta " + out_fasta,
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
	@summary: Change fasta headers to be compatible with picrust2. 
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
	@summary: Returns the excluded sequence, not insert into reference tree.
	@param fasta_file: [str] Path to the fasta file to process.
	@param tree_file: [str] Path to the tree file to process.
	@output: The file of no aligned sequence names.
	"""
	file = open(tree_file, "r")
	line = file.readline()
	list_cluster = re.findall("(Cluster_[0-9]+)", line)
	file.close()
	FH_input = FastaIO(in_fasta)
	excluded = open(excluded, "wt")

	no_excluded = True 
	for record in FH_input:
		if record.id not in list_cluster:
			excluded.write(record.id+"\n")
			no_excluded = False
	FH_input.close()
	if no_excluded:
		excluded.write('#No excluded OTUs.\n')
	excluded.close()

def restricted_float(in_arg):
	"""
	@summary: Custom argparse type to force an input float to be between 0 and 1.
	"""
	try:
		in_arg = float(in_arg)
	except ValueError:
		raise argparse.ArgumentTypeError(in_arg + " is not a floating-point "
										 "literal (i.e. not a proportion)")

	if in_arg < 0.0 or in_arg > 1.0:
		raise argparse.ArgumentTypeError(in_arg + "is not in range 0.0 - 1.0")
	return in_arg

def write_summary(in_fasta, align_out, biomfile, closest_ref_file, category, summary_file):
	"""
	@summary: Writes the process summary in one html file.
	@param summary_file: [str] path to the output html file.
	@param align_out: [str] path to the fasta file of unaligned OTU
	@param biomfile: [str] path to the input BIOM file.
	@param closest_ref_files: [str] Path to tmp colest ref file.
	@param category: ITS or 16S
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
	details_categorys =["FROGS Taxonomy","Picrust2 closest ID","Picrust2 closest reference name","Picrust2 closest taxonomy","Picrust2 closest distance from cluster","Comment", "Picrust2 closest reference sequence"]
	infos_otus = list()
	biom=BiomIO.from_json(biomfile)
	list_otu_all = []
	# record nb OTU and abundance
	for otu in FastaIO(in_fasta):
		list_otu_all.append(otu.id)
		number_otu_all +=1
		number_abundance_all += biom.get_observation_count(otu.id)

	if category == "16S":
		START_IMG_LINK = "<a href='https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
	elif category == "ITS":
		START_IMG_LINK = "<a href='https://mycocosm.jgi.doe.gov/"


	closest_ref = open(closest_ref_file)
	for li in closest_ref:
		li = li.strip().split('\t')
		if category in ['16S','ITS']:
			id_cur = li[2]
			li[2] = START_IMG_LINK + id_cur + "'>" + id_cur + '</a>'
		if "Closest_ref_name" not in li:
			infos_otus.append({
				'name': li[0],
				'data': list(map(str,li[1:]))
				})

	# record details about removed OTU
	if align_out is not None:
		for li in open(align_out).readlines():
			if re.match("(Cluster_[0-9]+)",li):
				summary_info['nb_removed'] +=1
				summary_info['abundance_removed'] += biom.get_observation_count(li.strip())
	
	summary_info['nb_kept'] = number_otu_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep1_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
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
	group_input.add_argument('-r', '--ref_dir', help='Directory containing reference sequence files, if studied marker is not 16S. (ex: $PICRUST2_PATH/default_files/fungi/fungi_ITS')
	group_input.add_argument('-p', '--placement_tool', default='epa-ng', help='Placement tool to use when placing sequences into reference tree. One of "epa-ng" or "sepp" must be input')
	group_input.add_argument('--min_align', type=restricted_float, default=0.8, help='Proportion of the total length of an input query ''sequence that must align with reference sequences. ''Any sequences with lengths below this value after ''making an alignment with reference sequences will ''be excluded from the placement and all subsequent ''steps. (default: %(default)d).')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--out_tree', default='FPStep1.tree', help='Tree output with insert sequences (format: newick).')
	group_output.add_argument('-e', '--excluded', default='FPSTep1_excluded.txt', help='List of sequences not inserted in the tree.')
	group_output.add_argument('-s', '--insert_fasta', default='FPStep1.fasta', help='sequences file without non insert sequences. (format: FASTA). [Default: %(default)s]')
	group_output.add_argument('-m', '--insert_biom', default='FPStep1.biom', help='abundance file without non insert sequences. (format: BIOM) [Default: %(default)s]')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--html', default='FPStep1_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.out_tree)[0])

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		if args.ref_dir is None or 'pro_ref' in args.ref_dir:
			category = '16S'
		else:
			category = 'ITS'

		Logger.static_write(args.log_file,'\n# Cleaning fasta headers\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
		tmp_fasta = tmp_files.add('cleaned.fasta')
		convert_fasta(args.input_fasta,tmp_fasta)

		tmp_place_seqs = tmp_files.add( 'tmp_place_seqs.log' )
		PlaceSeqs(tmp_fasta, args.out_tree, args.placement_tool, args.ref_dir, args.min_align, tmp_place_seqs).submit(args.log_file)
		# parse place_seqs.py output in order to retrieve references sequences alignment, necessary for find_closest_ref_sequences step.
		ref_aln = open(tmp_place_seqs).readlines()[0].strip().split()[4]


		tmp_multi_affiliations = tmp_files.add( 'multi_affiliations.tsv' )
		tmp_multi_affi_log = tmp_files.add( 'multi_affiliations.log' )
		multiAffiFromBiom(args.input_biom, tmp_multi_affiliations, tmp_multi_affi_log).submit(args.log_file)

		excluded_sequence(args.out_tree,args.input_fasta,args.excluded)

		RemoveSeqsBiomFasta(tmp_fasta, args.input_biom, args.insert_fasta, args.insert_biom, args.excluded).submit(args.log_file)

		closest_ref_files = tmp_files.add( "closest_ref.tsv" )
		tmp_find_closest_ref = tmp_files.add( 'tmp_find_closest_ref.log' )
		FindClosestsRefSequences(args.out_tree, args.input_biom, tmp_multi_affiliations, tmp_fasta, ref_aln, args.insert_biom, closest_ref_files).submit(args.log_file)
		write_summary(tmp_fasta, args.excluded, args.input_biom, closest_ref_files, category, args.html)

	finally:
		if not args.debug:
			tmp_files.deleteAll()
