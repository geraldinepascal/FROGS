#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

#Import
import os
import re
import sys
import json
import glob
import gzip
import shutil
import argparse


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

DESCRIPTION_DIR = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/description_mapfiles/")

##################################################################################################################################################
#
# COMMAND LINES 
#
##################################################################################################################################################
class PathwayPipeline(Cmd):
	"""
	@summary: pathway_pipeline.py : Infer the presence and abundances of pathways based on gene family abundances in a sample.
	"""
	def __init__(self, input_file, per_sequence_contrib, per_sequence_abun, per_sequence_function, pathways_abund, pathways_contrib, pathways_predictions, log):
		opt_contrib = ''
		if per_sequence_contrib:
			opt_contrib = ' --per_sequence_contrib --per_sequence_abun ' +  per_sequence_abun + ' --per_sequence_function ' + per_sequence_function

		Cmd.__init__(self,
				 'pathway_pipeline.py ',
				 'predict abundance pathway', 
				  " --input " + input_file + " --out_dir ./ " + opt_contrib + ' 2> ' + log,
				"--version")

		self.pathways_abund = pathways_abund
		self.per_sequence_contrib = per_sequence_contrib
		self.pathways_contrib = pathways_contrib
		self.pathways_predictions = pathways_predictions
		
	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip()

	def parser(self, log_file):
		with gzip.open('path_abun_unstrat.tsv.gz', 'rb') as f_in:
			with open(self.pathways_abund, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('path_abun_unstrat.tsv.gz')
		if self.per_sequence_contrib:
			with gzip.open('path_abun_contrib.tsv.gz', 'rb') as f_in:
				with open(self.pathways_contrib, 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
				os.remove('path_abun_contrib.tsv.gz')
			with gzip.open('path_abun_predictions.tsv.gz', 'rb') as f_in:
				with open(self.pathways_predictions, 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
				os.remove('path_abun_predictions.tsv.gz')

class AddDescriptions(Cmd):
	"""
	@summary: Adds a description column to a function abundance table and outputs a new file.
	"""
	def __init__(self, abund_file, description_file, out_file):
		"""
		@param function_file: [str] Path to input pathway file. (ex: pathway_out/path_abun_unstrat.tsv.gz)
		"""
		Cmd.__init__(self,
			'add_descriptions.py ',
			'Adds a description column to a function abundance table.',
			'--input ' + abund_file + ' --custom_map_table ' + description_file + ' --output ' + out_file,
			'--version')

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def formate_description_file(description_dir, description_out):
	"""
	@summary: concatane all picrust2 descriptions file into one temporary global description file.
	"""
	with open(description_out, 'wb') as outfile:
		for filename in glob.glob(description_dir+'/*'):
			if filename != description_out and not 'KEGG' in filename:
				with open(filename, 'rb') as readfile:
					shutil.copyfileobj(readfile, outfile)

def write_summary(strat_file, summary_file):
	"""
	@summary: Writes the process summary in one html file.
	@param summary_file: [str] path to the output html file.
	@param align_out: [str] path to the fasta file of unaligned OTU
	@param biomfile: [str] path to the input BIOM file.
	@param closest_ref_files: [str] Path to tmp colest ref file.
	@param category: ITS or 16S
	"""
	# to summary OTUs number && abundances number			   
	infos_otus = list()
	details_categorys =["Pathway", "Description" ,"Observation_sum"]
	START_METACYC_LINK = "<a href='https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object="

	abund = open(strat_file)
	for li in abund:
		if "pathway" in li:
			li = li.strip().split('\t')
			for sample in li[3:]:
				details_categorys.append(sample)
			break

	for li in abund:
		li = li.strip().split('\t')
		pathway = li[0]
		li[0] = START_METACYC_LINK + pathway + "'>" + pathway + '</a>'

		infos_otus.append({
			'name': li[0],
			'data': list(map(str,li[1:]))
			})
	# record details about removed OTU

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep4_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
		if "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_otus) )
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
	parser.add_argument('--skip_descriptions', default=False, action='store_true', help='Skipping add_descriptions.py step that add a description column to the function abundance table. (default: False')
	parser.add_argument('--per_sequence_contrib', default=False, action='store_true', help='Flag to specify that MinPath is run on the genes contributed by each sequence individualy. (in contrast to the default stratified output, which is the contribution to the community-wide pathway abundances.) Options --per_sequence_abun and --per_sequence_function need to be set when this option is used (default: False) ')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-i', '--input_file', required=True, type=str, help='Input TSV table of gene family abundances (either ''the stratified output of ' 'metagenome_pipeline.py).')
	group_input.add_argument('--per_sequence_abun', default=None, help='Path to table of sequence abundances across samples normalized by marker copy number (typically the normalized sequence abundance table output at the metagenome pipeline step: seqtab_norm.tsv by default). This input is required when the --per_sequence_contrib option is set. (default: None).')
	group_input.add_argument('--per_sequence_function', default=None, help='Path to table of function abundances per sequence, which was outputted at the hidden-state prediction step. This input is required when the --per_sequence_contrib option is set. Note that this file should be the same input table as used for the metagenome pipeline step (default: None).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-o', '--pathways_abund', default='FPStep4_path_abun_unstrat.tsv', help='Pathway abundance file output.')
	group_output.add_argument('--pathways_contrib', default=None, help='Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome.')
	group_output.add_argument('--pathways_predictions', default=None, help='Stratified output corresponding to contribution of predicted gene family abundances within each predicted genome.')
	group_output.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
	group_output.add_argument('-t', '--html', default='FPStep4_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )	
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.input_file)[0])
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		if args.per_sequence_contrib:
			if args.per_sequence_abun == None or args.per_sequence_function == None:
				parser.error("\n\n#ERROR : --per_sequence_abun and --per_sequence_function required when --per_sequence_contrib option is set!\n\n")
			if args.pathways_contrib is None:
				args.pathways_contrib = 'FPStep4_path_abun_contrib.tsv'
			if args.pathways_predictions is None:
				args.pathways_predictions = 'FPSTep4_path_abun_predictions.tsv'


		if (args.per_sequence_abun is not None or args.per_sequence_function is not None) and not args.per_sequence_contrib:
			parser.error("\n\n#ERROR : --per_sequence_contrib required when --per_sequence_contrib and --per_sequence_function option is set!\n\n")

		tmp_pathway = tmp_files.add( 'pathway_pipeline.log' )
		PathwayPipeline(args.input_file, args.per_sequence_contrib, args.per_sequence_abun, args.per_sequence_function, args.pathways_abund, args.pathways_contrib, args.pathways_predictions, tmp_pathway).submit(args.log_file)

		if not args.skip_descriptions:
			tmp_description_file = tmp_files.add('descriptions_file.tsv.gz')
			formate_description_file(DESCRIPTION_DIR, tmp_description_file )

			AddDescriptions(args.pathways_abund,  tmp_description_file, args.pathways_abund).submit( args.log_file)

		write_summary(args.pathways_abund, args.html)

	finally:
		if not args.debug:
			tmp_files.deleteAll()

