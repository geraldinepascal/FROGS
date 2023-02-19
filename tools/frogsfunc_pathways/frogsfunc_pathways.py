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

__author__ = ' Moussa Samb & Vincent Darbot & Geraldine Pascal - GENPHYSE '
__copyright__ = 'Copyright (C) 2022 INRAE'
__license__ = 'GNU General Public License'
__version__ = '4.0.1'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'dev'

import os
import sys
import json
import argparse
import pandas as pd

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
if os.getenv('DESCRIPTION_FILE'):
   DESCRIPTION_FILE=os.environ['DESCRIPTION_FILE']  
else:
   DESCRIPTION_FILE=os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "frogsfunc_suppdata/pathways_description_file.txt.gz"))
if os.getenv('PATHWAYS_HIERARCHY_FILE'):
	
   PATHWAYS_HIERARCHY_FILE =os.environ['PATHWAYS_HIERARCHY_FILE']  
else:
   PATHWAYS_HIERARCHY_FILE =os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "frogsfunc_suppdata/pathways_hierarchy.tsv"))

#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES 
#
##################################################################################################################################################
class PathwayPipeline(Cmd):
	"""
	@summary: pathway_pipeline.py : Infer the presence and abundances of pathways based on gene family abundances in a sample.
	"""
	def __init__(self, input_file, map_file, per_sequence_contrib, per_sequence_abun, per_sequence_function, output_dir, log):
		"""
		@param input_file: [str] Input TSV table of gene family abundances (frogsfunc_genefamilies_pred_metagenome_unstrat.tsv from frogsfunc_genefamilies.py.
		@param map_file: [str] Mapping file of pathways to reactions, necessary if marker studied is not 16S.
		@param per_sequence_contrib: [boolean] Flag to specify that MinPath is run on the genes contributed by each sequence individualy.
		@param per_sequence_abun: [str] Path to table of sequence abundances across samples normalized by marker copy number (if per_sequence_contrib).
		@param per_sequence_function: [str] Path to table of function abundances per sequence, which was outputted at the hidden-state prediction step (if per_sequence_contrib).
		@param pathways_abund: [str] Pathway abundance file output..
		@param pathways_contrib: [str] Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome (if per_sequence_contrib).
		@param pathways_predictions: [str] Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome.
		@param pathways_abund_per_seq: [str] Pathway abundance file output per sequences (if --per_sequence_contrib set).
		"""	
		opt = ''
		if per_sequence_contrib:
			opt = ' --per_sequence_contrib --per_sequence_abun ' +  per_sequence_abun + ' --per_sequence_function ' + per_sequence_function 
		if map_file is not None:
			opt += " --map " + map_file
			if os.path.basename(map_file) == "KEGG_pathways_to_KO.tsv" :
				opt += " --no_regroup "

		Cmd.__init__(self,
				 'pathway_pipeline.py ',
				 'predict abundance pathway', 
				  " --input " + input_file + " --out_dir " + output_dir + opt + ' 2> ' + log,
				"--version")
		
	def get_version(self):
		 return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

class ParsePathwayPipeline(Cmd):
	"""
	@summary: Parse results of PICRUSt2 pathway_pipeline.py software to rerieve additional informations (i.g. databases functions links)
	"""
	def __init__(self, out_dir, out_abund, per_sequence_contrib, contrib, predictions, abund_per_seq, log):
		opt = ''
		if per_sequence_contrib:
			opt += " --per-sequence-contrib --output-contrib " + contrib + " --output-predictions " + predictions + " --output-abund-per-seq " + abund_per_seq 	
		Cmd.__init__( self,
					  'frogsFuncUtils.py',
					  'Parse pathway_pipeline.py outputs.',
					  "parse-pathway --input-dir " + out_dir + " --output-abund " + out_abund + opt + " 2>> " + log,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

class Biom2tsv(Cmd):
	"""
	@summary: Converts BIOM file to TSV file.
	"""
	def __init__(self, in_biom, out_tsv):

		Cmd.__init__( self,
					  'biom2tsv.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-file " + in_biom + " --output-file " + out_tsv + " --fields @observation_name @sample_count" ,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

class Tsv2biom(Cmd):
	"""
	@summary: In order to creates a temporary biom file that links every pathway to samples abundances.
	This is necessary in order to display sunburst plots.

	"""
	def __init__(self, in_tsv, out_biom):

		Cmd.__init__( self,
					  'tsv_to_biom.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-tsv " + in_tsv + " --output-biom " + out_biom,
					  '--version' )

		self.in_tsv = in_tsv

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip() 

	def parser(self, log_file):
		f_in = pd.read_csv(self.in_tsv, sep='\t')
		sum_col = f_in.pop("observation_sum")
		f_in.to_csv(self.in_tsv ,sep='\t' ,index=False)

class FormateAbundances(Cmd):
	"""
	@summary: Formate pathway abundances file in order to add function classifications and display sunbursts graphs.
	"""
	def __init__(self, in_abund, tmp_sunburst, tmp_unstrat, hierarchy_file, log):

		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Formate pathway abundances file.',
			'formate-abundances --input-abundances ' + in_abund + ' --input-tmp-sunburst ' + tmp_sunburst + ' --input-tmp-unstrat ' + tmp_unstrat + ' --hierarchy-file ' + hierarchy_file + ' 2>> ' + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class TaxonomyTree(Cmd):
	"""
	@summary: Produces a tree with pathway abundances by sample in extended newick format.
	"""
	def __init__(self, in_biom, taxonomy_tag, out_tree, out_ids):
		"""
		@param in_biom: [str] The processed BIOM path.
		@param taxonomy_tag: [str] The metadata title for the taxonomy in BIOM file.
		@param out_tree: [str] Path to the enewick output.
		@param out_ids: [str] Path to the IDs/samples output.
		"""
		# Cmd
		Cmd.__init__( self,
					  'biomTools.py',
					  'Produces a taxonomy tree with counts by sample.',
					  'treeCount --input-file ' + in_biom + ' --taxonomy-key "' + taxonomy_tag + '" --output-enewick ' + out_tree + ' --output-samples ' + out_ids,
					  '--version' )
					  
	def get_version(self):   
		return Cmd.get_version(self, 'stdout').strip()  

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################


def check_basename_files(arg_name, file_path):
	'''
	Test if output file name specified by the user only contains the file name, without directory.
	'''
	if not os.path.basename(file_path) == file_path:
		return raise_exception( Exception( "\n\n#ERROR : --" + arg_name.replace('_','-') + \
		" should only contain a filename, without directory (You specified " + arg_value +" ). Please use --output-dir to specify the output directory.\n\n"))

def formate_input_file(input_file, tmp_tsv):
	"""
	@summary: formate gene abundances file in order to use pathways_pipeline.py
	"""
	FH_in = open(input_file).readlines()
	FH_out = open(tmp_tsv, 'wt')
	header = FH_in[0].strip().split('\t')
	formatted_col = ['classification', 'observation_sum', 'db_link']
	to_keep = list()
	to_write = list()
	for col_n in range(len(header)):
		if header[col_n] == "observation_name":
			header[col_n] = "function"
		if header[col_n] not in formatted_col:
			to_write.append(header[col_n])
			to_keep.append(col_n)
	FH_out.write('\t'.join(to_write)+'\n')
	for li in FH_in[1:]:
		li = li.strip().split('\t')
		to_write = list()
		for i in range(len(li)):
			if i in to_keep:
				to_write.append(li[i])
		FH_out.write('\t'.join(to_write)+'\n')

def normalized_abundances_file( strat_file ):
	"""
	@summary normalized final output table in order to make the values comparable between samples. 
	Normalisation is done as follow :
	 (value)/(sum of the value in that column)*10^6, which gives CPM values
	 @param strat_file: [str] path to frogsfunc_pathways output abundances file.
	"""
	df = pd.read_csv(strat_file,sep='\t')
	for column in df.iloc[:,3:]:
		df[column] = df[column] / df[column].sum() * 1000000
		df[column] = df[column].round(0).astype(int)
	df.to_csv(strat_file, sep='\t', index=False)

def write_summary(strat_file, tree_count_file, tree_ids_file, summary_file):
	"""
	@summary: Writes the process summary in one html file.
	@param tree_count_file [str]: newick file of pathway abudances per samples and per hierarchy.
	@param tree_ids_file: [str] file that link id to its sample.
	@param summary_file: [str] path to the output html file.
	"""
	# Get taxonomy distribution
	FH_tree_count = open( tree_count_file )
	newick_tree = FH_tree_count.readline()
	FH_tree_count.close()
	ordered_samples_names = list()
	FH_tree_ids = open( tree_ids_file )
	for line in FH_tree_ids:
		id, sample_name = line.strip().split( "\t", 1 )
		ordered_samples_names.append( sample_name )
	FH_tree_ids.close()

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "frogsfunc_pathways_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###TAXONOMIC_RANKS###" in line:
			line = line.replace( "###TAXONOMIC_RANKS###", json.dumps(HIERARCHY_RANKS) )
		elif "###SAMPLES_NAMES###" in line:
			line = line.replace( "###SAMPLES_NAMES###", json.dumps(ordered_samples_names) )
		elif "###DATA_SAMPLE###" in line:
			line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
		elif "###TREE_DISTRIBUTION###" in line:
			line = line.replace( "###TREE_DISTRIBUTION###", json.dumps(newick_tree) )
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
	parser = argparse.ArgumentParser( description='Infer the presence and abundances of pathways based on gene family abundances in a sample.' )
	parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	parser.add_argument('--per-sequence-contrib', default=False, action='store_true', help='If stratified option is activated, a new table is built. It will contain the abundances of each function of each OTU in each sample. (in contrast to the default stratified output, which is the contribution to the community-wide pathway abundances.) Options --per-sequence-abun and --per-sequence-function need to be set when this option is used (default: False) ')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-i', '--input-file', required=True, type=str, help='Input TSV function abundances table from FROGSFUNC_step3_function (unstratified table : frogsfunc_functions_unstrat.tsv).')
	group_input.add_argument('-m', '--map', type=str, help='File required if you are not analyzing 16S sequences with the Metacyc ("EC" function in the previous step) database. IF MARKER STUDYED STILL 16S: it must indicate the path to the PICRUSt2 KEGG pathways mapfile, if you chose "KO" in the previous step (the mapfile is available here : $PICRUSt2_PATH/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv) IF MARKER STUDYED IS ITS OR 18S: Path to mapping file of pathways to fungi reactions (the mapfile is available here : $PICRUSt2_PATH/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt ).')
	group_input.add_argument('--per-sequence-abun', default=None, help='Path to table of sequence abundances across samples normalized by marker copy number (typically the normalized sequence abundance table output at the metagenome pipeline step: frogsfunc_functions_marker_norm.tsv by default). This input is required when the --per-sequence-contrib option is set. (default: None).')
	group_input.add_argument('--per-sequence-function', default=None, help='Path to table of function abundances per sequence, which was outputted at the hidden-state prediction step (frogsfunc_copynumbers_predicted_functions.tsv by default). This input is required when the --per-sequence-contrib option is set. Note that this file should be the same input table as used for the metagenome pipeline step (default: None).')
	group_input.add_argument('--hierarchy-ranks', nargs='*', default=["Level1", "Level2", "Level3", "Pathway"], help='The ordered ranks levels used in the metadata hierarchy pathways. [Default: %(default)s]' )
	group_input.add_argument( '--normalisation', default=False, action='store_true', help='To normalise data after analysis. Values are divided by sum of columns , then multiplied by 10^6 (CPM values). [Default: %(default)s]')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-d', '--output-dir', default='frogsfunc_pathway_results', help='Output directory for pathway predictions.')
	group_output.add_argument('-o', '--output-pathways-abund', default='frogsfunc_pathways_unstrat.tsv', help='Pathway abundance file output. Default: %(default)s]')
	group_output.add_argument('--output-pathways-contrib', default=None, help='Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome.')
	group_output.add_argument('--output-pathways-predictions', default=None, help='Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome.')
	group_output.add_argument('--output-pathways-abund-per-seq', default=None, help='Pathway abundance file output per sequences (if --per-sequence-contrib set)')
	group_output.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
	group_output.add_argument('-t', '--summary', default='frogsfunc_pathways_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )	
	args = parser.parse_args()
	prevent_shell_injections(args)

	args_dict = vars(args)
	for arg_name, arg_value in args_dict.items():
		if arg_name.startswith('output') and arg_name != "output_dir" and arg_value is not None:
			check_basename_files(arg_name, arg_value)

	if args.per_sequence_contrib:
		if args.per_sequence_abun == None or args.per_sequence_function == None:
			parser.error("\n\n#ERROR : --per-sequence-abun and --per-sequence-function required when --per-sequence-contrib option is set!\n\n")
		if args.output_pathways_contrib is None:
			args.output_pathways_contrib = args.output_dir + '/frogsfunc_pathways_strat.tsv'
		if args.output_pathways_predictions is None:
			args.output_pathways_predictions = args.output_dir + '/frogsfunc_pathways_predictions.tsv'
		if args.output_pathways_abund_per_seq is None:
			args.output_pathways_abund_per_seq = args.output_dir + '/frogsfunc_pathways_unstrat_per_seq.tsv'

	if (args.per_sequence_abun is not None or args.per_sequence_function is not None) and not args.per_sequence_contrib:
		parser.error("\n\n#ERROR : --per-sequence-contrib required when --per-sequence-contrib and --per-sequence-function option is set!\n\n")

	args.output_pathways_abund = args.output_dir + "/" + args.output_pathways_abund
	args.summary = args.output_dir + "/" + args.summary
	tmp_files=TmpFiles(os.path.split(args.summary)[0])
	tmp_files_picrust =  TmpFiles(os.path.dirname(args.output_pathways_abund), prefix="")
	
	HIERARCHY_RANKS = ['Level1','Level2','Level3','Pathway']
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		tmp_pathway = tmp_files.add( 'pathway_pipeline.log' )
		tmp_tsv = tmp_files.add( 'genes_abundances_formatted.tsv')
		formate_input_file(args.input_file, tmp_tsv)
		##
		tmp_seqtab = tmp_files_picrust.add('path_abun_unstrat.tsv.gz')
		if args.per_sequence_contrib:
			tmp_contrib = tmp_files_picrust.add('path_abun_contrib.tsv.gz')
			tmp_predictions = tmp_files_picrust.add('path_abun_predictions.tsv.gz')
			tmp_unstrat_per_seq = tmp_files_picrust.add('path_abun_unstrat_per_seq.tsv.gz')
		##
		try:
			PathwayPipeline(tmp_tsv, args.map, args.per_sequence_contrib, args.per_sequence_abun, args.per_sequence_function, args.output_dir, tmp_pathway).submit(args.log_file)
		except:
			raise_exception( Exception("\n\n#Note that the default pathway and regroup mapfiles are meant for EC numbers with 16S sequences. KEGG pathways are not supported since KEGG is a closed-source database, but you can input custom pathway mapfiles with the flag --map, associated with the file available here: $PICRUSt2_PATH/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv. For ITS or 18S please use --map with the file available here: $PICRUSt2_PATH/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt. \n\n"))
			
		tmp_parse_pathway = tmp_files.add( 'parse_pathway.log' )

		ParsePathwayPipeline(args.output_dir, args.output_pathways_abund, args.per_sequence_contrib, args.output_pathways_contrib, args.output_pathways_predictions, args.output_pathways_abund_per_seq, tmp_parse_pathway).submit( args.log_file)

		tmp_formate_abundances = tmp_files.add( 'tmp_formate_abundances.log' )
		tmp_pathway_sunburst = tmp_files.add( "functions_unstrat_sunburst.tmp")
		tmp_pathway_unstrat = tmp_files.add( "functions_unstrat.tmp")
		FormateAbundances(args.output_pathways_abund, tmp_pathway_sunburst, tmp_pathway_unstrat, PATHWAYS_HIERARCHY_FILE, tmp_formate_abundances).submit( args.log_file)
		if args.normalisation:
			normalized_abundances_file( args.output_pathways_abund)
		tmp_biom = tmp_files.add( 'pathway_abundances.biom' )
		Tsv2biom( tmp_pathway_sunburst, tmp_biom ).submit( args.log_file)
		tree_count_file = tmp_files.add( "pathwayCount.enewick" )
		tree_ids_file = tmp_files.add( "pathwayCount_ids.tsv" )
		hierarchy_tag = "classification"
		TaxonomyTree( tmp_biom, hierarchy_tag, tree_count_file, tree_ids_file ).submit( args.log_file )

		write_summary( args.output_pathways_abund, tree_count_file, tree_ids_file, args.summary )

	finally:
		if not args.debug:
			tmp_files.deleteAll()
			tmp_files_picrust.deleteAll()

