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
import glob
import json
import shutil
import argparse
from tempfile import gettempdir

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class MetagenomePipeline(Cmd):
	"""
	@summary: Per-sample metagenome functional profiles are generated based on the predicted functions for each study sequence.
	"""
	def __init__(self, in_biom, marker, function, max_nsti, min_reads, min_samples, strat_out, function_abund, seqtab, weighted, contrib, log):
		"""
		@param: in_biom: Path to BIOM input file used in FPStep1
		"""
		opt = ' --strat_out ' if strat_out else ''

		Cmd.__init__(self,
				 'metagenome_pipeline.py ',
				 'Per-sample functional profiles prediction.', 
				 " --input " +  in_biom + " --marker " + marker + " --function " + function + " --out_dir ./ --max_nsti " + str(max_nsti) + " --min_reads " + str(min_reads) + " --min_samples " + str(min_samples) + opt + ' 2> ' + log,
				"--version")

		self.abund = function_abund
		self.seqtab = seqtab
		self.weighted = weighted
		self.strat = strat_out
		self.contrib = contrib

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

	def parser(self, log_file):
		with gzip.open('pred_metagenome_unstrat.tsv.gz', 'rb') as f_in:
			with open(self.abund, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('pred_metagenome_unstrat.tsv.gz')
		with gzip.open('seqtab_norm.tsv.gz', 'rb') as f_in:
			with open(self.seqtab, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('seqtab_norm.tsv.gz')
		with gzip.open('weighted_nsti.tsv.gz', 'rb') as f_in:
			with open(self.weighted, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('weighted_nsti.tsv.gz')
		if self.strat:
			with gzip.open('pred_metagenome_contrib.tsv.gz', 'rb') as f_in:
				with open(self.contrib, 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
				os.remove('pred_metagenome_contrib.tsv.gz')

class AddDescriptions(Cmd):
	"""
	@summary: Adds a description column to a function abundance table and outputs a new file.
	"""
	def __init__(self, function_file, description_file, out_file, log):
		"""
		@param function_file: [str] Path to input function abundance table. (ex: EC_metagenome_out/pred_metagenome_unstrat.tsv.gz)
		"""
		Cmd.__init__(self,
			'add_descriptions.py ',
			'Adds a description column to a function abundance table.',
			'--input ' + function_file + ' --custom_map_table ' + description_file + ' --output ' + out_file + " 2> " + log,
			'--version')

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

class Biom2tsv(Cmd):
	"""
	@summary: Converts BIOM file to TSV file.
	"""
	def __init__(self, in_biom, out_tsv):

		Cmd.__init__( self,
					  'biom2tsv.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-file " + in_biom + " --output-file " + out_tsv,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip() 

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

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
	details_categorys =["Function", "Description" ,"Observation_sum"]
	START_GENBANK_LINK = "<a href='https://www.genome.jp/dbget-bin/www_bget?"
	START_COG_LINK = "<a href='https://www.ncbi.nlm.nih.gov/research/cog/cog/"
	START_PFAM_LINK = "<a href='https://pfam.xfam.org/family/"
	START_TIGR_LINK = "<a href='https://0-www-ncbi-nlm-nih-gov.linyanti.ub.bw/genome/annotation_prok/evidence/"

	abund = open(strat_file)
	for li in abund:
		if "function" in li:
			li = li.strip().split('\t')
			for sample in li[3:]:
				details_categorys.append(sample)
			break

	for li in abund:
		li = li.strip().split('\t')
		function = li[0]
		if "COG" in function:
			li[0] = START_COG_LINK + function + "'>" + function + '</a>'
		if "PF" in function:
			li[0] = START_PFAM_LINK + function + "'>" + function + '</a>'
		if "TIGR" in function:
			li[0] = START_TIGR_LINK + function + "'>" + function + '</a>'
		elif re.search('K[0-9]{5}',function) or "EC:" in function:
			li[0] = START_GENBANK_LINK + function + "'>" + function + '</a>'

		for i in range(len(li[2:])):
			li[i+2] = round(float(li[i+2]),1)

		infos_otus.append({
			'name': li[0],
			'data': list(map(str,li[1:]))
			})
	# record details about removed OTU

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep3_tpl.html") )
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
	parser = argparse.ArgumentParser( description='Per-sample functional profiles prediction.' )
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	parser.add_argument('--strat_out', default=False, action='store_true', help='Output table stratified by sequences as well. By default this will be in \"contributional\" format ''(i.e. long-format) unless the \"--wide_table\" ''option is set. The startified outfile is named ''\"metagenome_contrib.tsv.gz\" when in long-format.')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input_biom', required=True, type=str, help='Input table of sequence abundances (BIOM file used in FPStep1).')
	group_input.add_argument('-f', '--function', required=True, type=str, help='Table of predicted gene family copy numbers ''(FPStep2 output, ex FPStep2_all_predicted.tsv).')
	group_input.add_argument('-m', '--marker', required=True, type=str, help='Table of predicted marker gene copy numbers ''(FPStep2 output, ex FPStep2_marker_nsti_predicted.tsv.')
	group_input.add_argument('--add_description', default=False, action='store_true', help='Flag to adds a description column to the function abundance table')
	group_input.add_argument('--max_nsti', type=float, default=2.0, help='Sequences with NSTI values above this value will ' 'be excluded (default: %(default)d).')
	group_input.add_argument('--min_reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for ''each input ASV. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	group_input.add_argument('--min_samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an ASV needs to be ''identfied within. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('--function_abund', default='FPStep3_pred_metagenome_unstrat.tsv', help='Output file for metagenome predictions abundance. (default: %(default)s).')
	group_output.add_argument('--seqtab', default='FPStep3_seqtab_norm.tsv', help='This output file will contain abundance normalized. (default: %(default)s).')
	group_output.add_argument('--weighted', default='FPStep3_weighted_nsti.tsv', help='This output file will contain the nsti value per sample (format: TSV). [Default: %(default)s]' )
	group_output.add_argument('--contrib', default=None, help=' Stratified output that represents contributions to community-wide abundances (ex pred_metagenome_contrib.tsv)')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--html', default='FPStep3_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )	
	args = parser.parse_args()
	prevent_shell_injections(args)

	if args.strat_out and args.contrib is None:
		args.contrib = 'FPStep3_pred_metagenome_contrib.tsv'

	if not args.strat_out and args.contrib is not None:
		parser.error('--contrib FILENAME must be include with --strat_out flag')

	tmp_files=TmpFiles(os.path.split(args.marker)[0])
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		#temp tsv file necessary for metagenome_pipeline.py
		tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
		Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )

		tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )	
		MetagenomePipeline(tmp_biom_to_tsv, args.marker, args.function, args.max_nsti, args.min_reads, args.min_samples, args.strat_out, args.function_abund, args.seqtab, args.weighted, args.contrib, tmp_metag_pipeline).submit( args.log_file )

		if args.add_description != None:
			description_file = 'default_files/pathways_description_file.txt.gz'
			tmp_add_descriptions = tmp_files.add( 'tmp_add_descriptions.log' )	
			pred_file = args.function_abund
			AddDescriptions(pred_file,  description_file, pred_file, tmp_add_descriptions).submit( args.log_file)

		write_summary(args.function_abund, args.html)
	finally:
		if not args.debug:
			tmp_files.deleteAll()