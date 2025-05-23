#!/usr/bin/env python3

__author__ = 'Moussa Samb - GENPHYSE & Vincent Darbot - GENPHYSE & Geraldine Pascal - GENPHYSE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'prod'

import os
import re
import sys
import json
import math
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
	@summary: place cluster sequences (i.e. ASVs) into a reference tree.
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
		'Place studies sequences (i.e. ASVs) on reference tree.',
		'--study_fasta ' + in_fasta + ' --out_tree ' + out_tree + ' --placement_tool ' + placement_tool + " --min_align " + str(min_align) + opt + " --verbose 2>> " + log,
		'--version')

	def get_version(self):
		return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

class FindClosestsRefSequences(Cmd):
	'''
	@summary: find ASVs closest reference sequences into a reference tree. 
	'''
	def __init__(self, in_tree, in_biom, in_fasta, ref_aln, out_biom, out_summary, log):
		'''
		@param in_tree: [str] Path to resulting tree file with insert clusters sequences.(place_seqs.py output).
		@param in_biom: [str] Path to BIOM input file.
		@param in_fasta: [str]	 Path to input fasta file of unaligned cluster sequences.
		@param ref_aln [str]: Path to the alignment file of reference sequences.
		@param out_biom [str]: Path to output Biom file with PICRUSt2 taxonomic affiliations informations.
		@summary out_summary: [str] Path to output summary file.
		@note : Header of summary file:
		Cluster	FROGS Taxonomy	PICRUSt2_closest_ID	PICRUSt2_closest_reference_name	PICRUSt2_closest_taxonomy	
		PICRUSt2_closest_distance_from_cluster_(NSTI)	FROGS_and_PICRUSt2_lowest_same_taxonomic_rank	Comment	Cluster_sequence	PICRUSt2_closest_reference_sequence
		'''
		Cmd.__init__(self,
			'find_closest_ref_sequence.py',
			'find ASVs closests reference sequences into a reference tree.',
			'--input-tree ' + in_tree + ' --input-biom ' + in_biom + ' --input-fasta ' + in_fasta + ' --ref-aln ' + ref_aln + ' --output-biom ' + out_biom + ' --output-tsv ' + out_summary + " --log-file " + log,
			'--version')
		self.log_file = log

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

class HspMarker(Cmd):
	"""
	@summary: Predict number of marker copies (16S, 18S or ITS) for each cluster sequence (i.e ASV).
	"""
	def __init__(self, tree, marker_type, marker_table, hsp_method, biom_file, output, log, is_debug):
		"""
		@param observed_marker_table: [str] Path to marker table file if marker studied is not 16S.
		@param in_tree: [str] Path to resulting tree file with inserted clusters sequences from frogsfunc_placeseqs.
		@param hsp_method: [str] HSP method to use.
		@param output: [str] PICRUSt2 marker output file.
		"""
		debug = ""
		if is_debug:
			debug = " --debug "
		if marker_type != "16S":
			opt = ' --input-marker-table ' + marker_table
		else:
			opt = ''

		Cmd.__init__(self,
				 'launch_hsp.py',
				 'predict marker copy number per sequence.', 
				  debug + ' marker --input-tree ' + tree + ' --marker-type ' + marker_type + opt + ' --hsp-method ' + hsp_method + ' -o ' + output + '  2> ' + log,
				"--version")

		self.output = output
		self.biom_file = biom_file

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()
	
	def parser(self, log_file):
		biom = BiomIO.from_json(self.biom_file)
		with open(self.output) as fi:
			for li in fi:
				if "metadata_NSTI" in li:
					continue
				li = li.strip().split('\t')
				cluster = li[0]
				NSTI = li[2]
				biom.add_metadata(cluster, "NSTI", NSTI, "observation", erase_warning = False)
		BiomIO.write(self.biom_file, biom)

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

def rounding(nb):
	'''
	@summary: Rounding numbers decimal 
	'''
	if re.search("^[0-9]{1}[.][0-9]+e",str(nb)):
		start = re.compile("[0-9][.][0-9]{1,2}")
		end = re.compile("e-[0-9]+")
		return float("".join(start.findall(str(nb))+end.findall(str(nb))))

	elif re.search("[0][.][0-9]+",str(nb)):
		return(round(nb,2))

	elif re.search("[0][.][0]+",str(nb)):
		motif = re.compile("[0][.][0]+[0-9]{2}")
		return float("".join(motif.findall(str(nb))))

	else:
		return(round(nb,2))

def write_summary(in_fasta, excluded_file, biomfile, closest_ref_file, category, log_find_closest_ref, depth_nsti_file, summary_file):
	"""
	@param in_fasta: [str] path to the input fasta file.
	@param excluded_file: [str] List of excluded sequences from the reference tree.
	@param biomfile: [str] path to the input BIOM file.
	@param closest_ref_files: [str] Path to closests reference information file (find_closest_ref_sequence.py output).
	@param category: ITS or 16S
	@param summary_file: [str] path to the output html file.
	"""
	# to summary ASVs number && abundances number			   
	summary_info = {
	   'nb_kept' : 0,
	   'nb_removed' : 0,
	   'abundance_kept' : 0,
	   'abundance_removed' : 0	   
	}
	number_asv_all = 0
	number_abundance_all = 0

	details_categorys =["Nb sequences","FROGS Taxonomy","PICRUSt2 closest ID (JGI)","PICRUSt2 closest reference name","PICRUSt2 closest taxonomy","NSTI", "NSTI Confidence" ,"Lowest same taxonomic rank between FROGS and PICRUSt2","Comment"]
	infos_asvs = list()
	biom=BiomIO.from_json(biomfile)
	list_asv_all = list()
	# record nb ASV and abundance
	for asv in FastaIO(in_fasta):
		list_asv_all.append(asv.id)
		number_asv_all +=1
		number_abundance_all += biom.get_observation_count(asv.id)

	if category == "16S":
		START_IMG_LINK = "<a href='https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid="
	elif category == "ITS":
		START_IMG_LINK = "<a href='https://mycocosm.jgi.doe.gov/"
	
	FH_in = open(log_find_closest_ref)
	for li in FH_in:
		if li.startswith('#Max NSTI'):
			max_nsti = float(li.strip().split()[-1])
	max_nsti = math.ceil( max_nsti * 50 + 1) 

	closest_ref = open(closest_ref_file)
	FH_log = Logger( depth_nsti_file )
	step_nsti = [i/50 for i in range(0, max_nsti)]
	cluster_kept = dict()

	for cur_nsti in step_nsti:
		cluster_kept[cur_nsti] = { 'Nb' : 0, 'Abundances' : 0 }
	for li in closest_ref:
		if li.startswith('#ASV'):
			continue
		li = li.strip().split('\t')
		for cur_nsti in step_nsti[::-1]:
			if float(li[6]) <= cur_nsti:
				cluster_kept[cur_nsti]['Nb']+=1
				cluster_kept[cur_nsti]['Abundances']+=int(li[1])
			else:
				break

		picrust_id_cur = li[3]
		li[3] = START_IMG_LINK + picrust_id_cur + "'target=\"_blank\">" + picrust_id_cur + '</a>'
		infos_asvs.append({
			'name': li[0],
			'data': list(li[1:-1])
			})
	closest_ref.close()
				
	clusters_size = list()
	abundances_size = list()
	for nsti,clusters in cluster_kept.items():
		clusters_size.append(clusters['Nb'])
		abundances_size.append(clusters['Abundances'])
		FH_log.write("\t".join([str(nsti), str(clusters['Nb']), str(clusters['Abundances']) ])+"\n")

	clusters_size = sorted(clusters_size)
	abundances_size = sorted(abundances_size)
	total_abundances = abundances_size[-1]
	proportions = [ rounding( i / total_abundances * 100 ) for i in abundances_size]
	# record details about removed ASV
	FH_excluded = open(excluded_file, 'rt').readlines()
	for li in FH_excluded:
		if not li.startswith('#No excluded ASV.'):
			cluster = li.strip()
			summary_info['nb_removed'] +=1
			summary_info['abundance_removed'] += biom.get_observation_count(cluster)

	summary_info['nb_kept'] = number_asv_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "frogsfunc_placeseqs_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
		elif "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_asvs) )
		elif "###REMOVE_DATA###" in line:
			line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
		elif "###CLUSTERS_SIZES###" in line:
			line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
		elif "###ABUNDANCES_SIZES###" in line:
			line = line.replace( "###ABUNDANCES_SIZES###", json.dumps( abundances_size) )
		elif "###STEP_NSTI###" in line:
			line = line.replace( "###STEP_NSTI###", json.dumps(step_nsti) )
		elif "###FROGS_VERSION###" in line:
			line = line.replace( "###FROGS_VERSION###", "\""+str(__version__)+"\"" )
		elif "###FROGS_TOOL###" in line:
			line = line.replace( "###FROGS_TOOL###", "\""+ os.path.basename(__file__)+"\"" )
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
	parser = argparse.ArgumentParser(description="place studies sequences (i.e. ASVs) into a reference tree.")
	parser.add_argument('--version', action='version', version=__version__)
	parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('--input-fasta', required=True, help="Input fasta file of unaligned studies sequences.")
	group_input.add_argument('--input-biom', required=True, help='Input biom file of unaligned studies sequences.')
	group_input.add_argument('--ref-dir', help='If marker studied is not 16S, this is the directory containing reference sequence files (for ITS, see: $PICRUST2_PATH/default_files/fungi/fungi_ITS')
	group_input.add_argument('--placement-tool', default='epa-ng', choices=["epa-ng", "sepp"], help='Tool to place sequences into reference tree. Note that epa-ng is more sensitiv but very memory and computing power intensive. Warning : sepp is not usable for ITS and 18S analysis [Default: %(default)s]')
	group_input.add_argument('--min-align', type=restricted_float, default=0.8, help='Proportion of the total length of an input query sequence that must align with reference sequences. Any sequences with lengths below this value after making an alignment with reference sequences will be excluded from the placement and all subsequent steps. [Default: %(default)s].')
	group_input.add_argument('--input-marker-table',help="The input marker table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. (ex $PICRUSt2_PATH/default_files/fungi/ITS_counts.txt.gz).")	
	group_input.add_argument('--hsp-method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use. mp: predict discrete traits using max parsimony. emp_prob: predict discrete traits based on empirical state probabilities across tips. subtree_average: predict continuous traits using subtree averaging. pic: predict continuous traits with phylogentic independent contrast. scp: reconstruct continuous traits using squared-change parsimony [Default: %(default)s].')   
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('--output-tree', default='frogsfunc_placeseqs_tree.nwk', help='Reference tree output with insert sequences (format: newick). [Default: %(default)s]')
	group_output.add_argument('--excluded', default='frogsfunc_placeseqs_excluded.txt', help='List of sequences not inserted in the tree. [Default: %(default)s]')
	group_output.add_argument('--output-fasta', default='frogsfunc_placeseqs.fasta', help='Fasta file without non insert sequences. (format: FASTA). [Default: %(default)s]')
	group_output.add_argument('--output-biom', default='frogsfunc_placeseqs.biom', help='Biom file without non insert sequences. (format: BIOM) [Default: %(default)s]')
	group_output.add_argument('--closests-ref', default='frogsfunc_placeseqs_closests_ref_sequences.txt', help='Informations about Clusters (i.e ASVs) and PICRUSt2 closest reference from cluster sequences (identifiants, taxonomies, phylogenetic distance from reference, nucleotidics sequences). [Default: %(default)s]')
	group_output.add_argument('--html', default='frogsfunc_placeseqs_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	group_output.add_argument('--output-marker', default="frogsfunc_marker.tsv", type=str, help='Output table of predicted marker gene copy numbers per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped. [Default: %(default)s]')	
	group_output.add_argument('--log-file', default=sys.stdout, help='List of commands executed. [Default: stdout]')
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

		if category == "ITS" and args.input_marker_table is None:
			raise_exception( Exception ("\n\n#ERROR : --input-marker-table required when studied marker is not 16S!\n\n"))

		tmp_hsp_marker = tmp_files.add( 'tmp_hsp_marker.log' )
		HspMarker(args.output_tree, category, args.input_marker_table, args.hsp_method, args.output_biom, args.output_marker, tmp_hsp_marker, args.debug).submit(args.log_file)

		tmp_find_closest_ref = tmp_files.add( 'tmp_find_closest_ref.log' )
		FindClosestsRefSequences(args.output_tree, args.output_biom, args.output_fasta, ref_aln, args.output_biom, args.closests_ref, tmp_find_closest_ref).submit(args.log_file)
		tmp_depth_nsti = tmp_files.add( 'depth_nsti.txt' )
		write_summary(tmp_fasta, args.excluded, args.input_biom, args.closests_ref, category, tmp_find_closest_ref, tmp_depth_nsti, args.html)

	finally:
		if not args.debug:
			tmp_files.deleteAll()
