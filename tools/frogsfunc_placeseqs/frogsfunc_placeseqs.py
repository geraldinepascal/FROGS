#!/usr/bin/env python3
#
# Copyright (C) 2022 INRAE
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

__author__ = ' Moussa Samb & Vincent Darbot & Geraldine Pascal GENPHYSE '
__copyright__ = 'Copyright (C) 2022 INRAE'
__license__ = 'GNU General Public License'
__version__ = '4.0.1'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'dev'

import os
import re
import sys
import json
import argparse

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
	@summary: place cluster sequences (i.e. OTUs) into a reference tree.
	"""
	def __init__(self, in_fasta, out_tree, placement_tool, ref_dir, min_align, log):
		"""
		@param in_fasta: [str] Path to input fasta file of unaligned cluster sequences.
		@param out_tree: [str] Path to output resulting tree file with insert clusters sequences.
		@param placement_tool: [str] Placement tool to use (epa-ng or sepp).
		@param ref_dir: [str] Directory containing reference sequence files.
		@param min_align: [float] Proportion of the total length of an input query sequence that must align with reference sequences.
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
		return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

class FindClosestsRefSequences(Cmd):
	'''
	@summary: find OTUs closest reference sequences into a reference tree. 
	'''
	def __init__(self, in_tree, in_biom, in_fasta, ref_aln, out_biom, out_summary, log):
		'''
		@param in_tree: [str] Path to resulting tree file with insert clusters sequences.(place_seqs.py output).
		@param in_biom: [str] Path to BIOM input file.
		@param in_fasta: [str] 	Path to input fasta file of unaligned cluster sequences.
		@param ref_aln [str]: Path to the alignment file of reference sequences.
		@param out_biom [str]: Path to output Biom file with PICRUSt2 taxonomic affiliations informations.
		@summary out_summary: [str] Path to output summary file.
		@note : Header of summary file:
		Cluster	FROGS Taxonomy	PICRUSt2_closest_ID	PICRUSt2_closest_reference_name	PICRUSt2_closest_taxonomy	
		PICRUSt2_closest_distance_from_cluster_(NSTI)	FROGS_and_PICRUSt2_lowest_same_taxonomic_rank	Comment	Cluster_sequence	PICRUSt2_closest_reference_sequence
		'''
		Cmd.__init__(self,
			'find_closest_ref_sequence.py',
			'find OTUs closests reference sequences into a reference tree.',
			'--input-tree ' + in_tree + ' --input-biom ' + in_biom + ' --input-fasta ' + in_fasta + ' --ref-aln ' + ref_aln + ' --output-biom ' + out_biom + ' --output-tsv ' + out_summary + " 2>> " + log,
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
			'--input-biom ' + in_biom + ' --input-fasta ' + in_fasta + ' --excluded-sequences ' + excluded_file + ' --output-biom ' + out_biom + " --output-fasta " + out_fasta,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class ConvertFasta(Cmd):
	'''
	@summary: Change fasta headers to be compatible with PICRUSt2.
	'''
	def __init__(self, in_fasta, out_fasta, log):
		'''
		@param input_fasta: [str] Path to fasta input file.
		@param output_fasta: [str] Path to fasta output file.
		'''
		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Change fasta headers to be compatible with PICRUSt2.',
			'convert-fasta --input-fasta ' + in_fasta + ' --output-fasta ' + out_fasta + ' 2>> ' + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class ExcludedSequences(Cmd):
	'''
	@summary: Returns the excluded sequence, not insert into reference tree.
	'''
	def __init__(self, in_tree, in_fasta, out_excluded, log):
		'''
		@param in_tree: [str] 'PICRUSt2 output tree with inserts sequences.
		@param in_fasta: [str] Path to input fasta file.
		@param out_excluded: [str] Path to excluded file with clusters not inserts in the reference tree.
		'''
		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Change fasta headers to be compatible with PICRUSt2.',
			'excluded-sequences --input-fasta ' + in_fasta + ' --input-tree ' + in_tree + ' --output-excluded ' + out_excluded + ' 2>> ' + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

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

def write_summary(in_fasta, excluded_file, biomfile, closest_ref_file, category, depth_nsti_file, summary_file):
	"""
	@param in_fasta: [str] path to the input fasta file.
	@param align_out: [str] path to the fasta file of unaligned OTU
	@param biomfile: [str] path to the input BIOM file.
	@param closest_ref_files: [str] Path to closests reference information file (find_closest_ref_sequence.py output).
	@param category: ITS or 16S
	@param summary_file: [str] path to the output html file.
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
	details_categorys =["Nb sequences","FROGS Taxonomy","PICRUSt2 closest ID (JGI)","PICRUSt2 closest reference name","PICRUSt2 closest taxonomy","NSTI", "NSTI Confidence" ,"Lowest same taxonomic rank between FROGS and PICRUSt2","Comment"]
	infos_otus = list()
	biom=BiomIO.from_json(biomfile)
	list_otu_all = list()
	# record nb OTU and abundance
	for otu in FastaIO(in_fasta):
		list_otu_all.append(otu.id)
		number_otu_all +=1
		number_abundance_all += biom.get_observation_count(otu.id)

	if category == "16S":
		START_IMG_LINK = "<a href='https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
	elif category == "ITS":
		START_IMG_LINK = "<a href='https://mycocosm.jgi.doe.gov/"


	closest_ref = open(closest_ref_file).readlines()
	for li in closest_ref[1:]:
		li = li.strip().split('\t')
		if category in ['16S','ITS']:
			try:
				id_cur = li[3]
				li[3] = START_IMG_LINK + id_cur + "'target=\"_blank\">" + id_cur + '</a>'
				infos_otus.append({
					'name': li[0],
					'data': list(li[1:-2])
					})
			except:
				continue

	FH_log = Logger( depth_nsti_file )
	step_nsti = [i/50 for i in range(0,101)] 
	cluster_kept = dict()
	for cur_nsti in step_nsti:
		cluster_kept[cur_nsti] = { 'Nb' : 0, 'Abundances' : 0 }
		for li in closest_ref[1:]:
			li = li.strip().split('\t')
			if float(li[6]) <= cur_nsti:
				cluster_kept[cur_nsti]['Nb']+=1
				cluster_kept[cur_nsti]['Abundances']+=int(li[1])
				
	clusters_size = list()
	abundances_size = list()
	nstis = list()
	for nsti,clusters in cluster_kept.items():
		clusters_size.append(clusters['Nb'])
		abundances_size.append(clusters['Abundances'])
		nstis.append(float(nsti))
		FH_log.write("\t".join([str(nsti), str(clusters['Nb']), str(clusters['Abundances']) ])+"\n")

	nstis = sorted(nstis)
	clusters_size = sorted(clusters_size)
	abundances_size = sorted(abundances_size)
	# record details about removed OTU
	FH_excluded = open(excluded_file, 'rt').readlines()
	for li in FH_excluded:
		cluster = li .strip()
		summary_info['nb_removed'] +=1
		summary_info['abundance_removed'] += biom.get_observation_count(cluster)

	summary_info['nb_kept'] = number_otu_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "frogsfunc_placeseqs_tpl_test.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
		elif "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_otus) )
		elif "###REMOVE_DATA###" in line:
			line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
		elif "###CLUSTERS_SIZES###" in line:
			line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
		elif "###ABUNDANCES_SIZES###" in line:
			line = line.replace( "###ABUNDANCES_SIZES###", json.dumps( abundances_size) )
		elif "###NSTI_THRESH###" in line:
			line = line.replace( "###NSTI_THRESH###", json.dumps(nstis) )
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
	group_input.add_argument('-i', '--input-fasta', required=True, help="Input fasta file of unaligned studies sequences.")
	group_input.add_argument('-b', '--input-biom', required=True, help='Input biom file of unaligned studies sequences.')
	group_input.add_argument('-r', '--ref-dir', help='If marker studied is not 16S, this is the directory containing reference sequence files (for ITS, see: $PICRUST2_PATH/default_files/fungi/fungi_ITS')
	group_input.add_argument('-p', '--placement-tool', default='epa-ng', choices=["epa-ng", "sepp"], help='Tool to place sequences into reference tree. Note that epa-ng is more sensitiv but very memory and computing power intensive. Warning : sepp is not usable for ITS and 18S analysis [Default: %(default)s]')
	group_input.add_argument('--min-align', type=restricted_float, default=0.8, help='Proportion of the total length of an input query sequence that must align with reference sequences. Any sequences with lengths below this value after making an alignment with reference sequences will be excluded from the placement and all subsequent steps. (default: %(default)s).')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--output-tree', default='frogsfunc_placeseqs_tree.nwk', help='Reference tree output with insert sequences (format: newick). [Default: %(default)s]')
	group_output.add_argument('-e', '--excluded', default='frogsfunc_placeseqs_excluded.txt', help='List of sequences not inserted in the tree. [Default: %(default)s]')
	group_output.add_argument('-s', '--output-fasta', default='frogsfunc_placeseqs.fasta', help='Fasta file without non insert sequences. (format: FASTA). [Default: %(default)s]')
	group_output.add_argument('-m', '--output-biom', default='frogsfunc_placeseqs.biom', help='Biom file without non insert sequences. (format: BIOM) [Default: %(default)s]')
	group_output.add_argument('-c', '--closests-ref', default='frogsfunc_placeseqs_closests_ref_sequences.txt', help='Informations about Clusters (i.e OTUs) and PICRUSt2 closest reference from cluster sequences (identifiants, taxonomies, phylogenetic distance from reference, nucleotidics sequences). [Default: %(default)s]')
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--summary', default='frogsfunc_placeseqs_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.output_tree)[0])

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		if args.ref_dir is None or 'pro_ref' in args.ref_dir:
			category = '16S'
		else:
			category = 'ITS'

		if args.placement_tool == "sepp" and category == "ITS":
			raise_exception( Exception ("\n\n#ERROR : You can't use sepp for ITS and 18S analysis.\n\n" ))

		Logger.static_write(args.log_file,'\n# Cleaning fasta headers\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
		tmp_fasta = tmp_files.add('cleaned.fasta')
		tmp_convert_fasta = tmp_files.add( 'tmp_convert_fasta.log' )
		ConvertFasta(args.input_fasta, tmp_fasta, tmp_convert_fasta).submit(args.log_file)

		tmp_place_seqs = tmp_files.add( 'tmp_place_seqs.log' )
		PlaceSeqs(tmp_fasta, args.output_tree, args.placement_tool, args.ref_dir, args.min_align, tmp_place_seqs).submit(args.log_file)
		# parse place_seqs.py output in order to retrieve references sequences alignment, necessary for find_closest_ref_sequences step.
		ref_aln = open(tmp_place_seqs).readlines()[0].strip().split()[4]

		tmp_excluded = tmp_files.add( 'tmp_excluded.log' )
		ExcludedSequences(args.output_tree, args.input_fasta, args.excluded, tmp_excluded).submit(args.log_file)

		RemoveSeqsBiomFasta(tmp_fasta, args.input_biom, args.output_fasta, args.output_biom, args.excluded).submit(args.log_file)

		tmp_find_closest_ref = tmp_files.add( 'tmp_find_closest_ref.log' )
		FindClosestsRefSequences(args.output_tree, args.output_biom, args.output_fasta, ref_aln, args.output_biom, args.closests_ref, tmp_find_closest_ref).submit(args.log_file)
		tmp_depth_nsti = tmp_files.add( 'depth_nsti.txt' )
		write_summary(tmp_fasta, args.excluded, args.input_biom, args.closests_ref, category, tmp_depth_nsti, args.summary)

	finally:
		if not args.debug:
			tmp_files.deleteAll()
