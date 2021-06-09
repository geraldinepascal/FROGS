#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
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
import argparse
from numpy import median
import shutil
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

DESCRIPTION_DIR = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/description_mapfiles/")
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class MetagenomePipeline(Cmd):
	"""
	@summary: Per-sample metagenome functional profiles are generated based on the predicted functions for each study sequence.
	"""
	def __init__(self, in_biom, marker, function, out_dir, max_nsti, min_reads, min_samples, strat_out, log):
		"""
		@param: in_biom: Path to BIOM input file used in FPStep1
		"""
		opt = ' --strat_out ' if strat_out else ''
	
		Cmd.__init__(self,
				 'metagenome_pipeline.py ',
				 'Per-sample functional profiles prediction.', 
				 " --input " +  in_biom + " --marker " + marker + " --function " + function + " --out_dir " + out_dir + " --max_nsti " + str(max_nsti) + " --min_reads " + str(min_reads) + " --min_samples " + str(min_samples) + opt + ' 2> ' + log,
				"--version") 
	  
	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

class AddDescriptions(Cmd):
	"""
	@summary: Adds a description column to a function abundance table and outputs a new file.
	"""
	def __init__(self, function_file, description_file, out_file):
		"""
		@param function_file: [str] Path to input function abundance table. (ex: EC_metagenome_out/pred_metagenome_unstrat.tsv.gz)
		@param function_type: [str] Function type between COG,EC,KO,PFAM and TIGRFAM.
		"""
		Cmd.__init__(self,
			'add_descriptions.py ',
			'Adds a description column to a function abundance table.',
			'--input ' + function_file + ' --custom_map_table ' + description_file + ' --output ' + out_file,
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

def formate_description_file(description_dir, description_out):
	"""
	@summary: concatane all picrust2 descriptions file into one temporary global description file.
	"""
	with open(description_out, 'wb') as outfile:
		for filename in glob.glob(description_dir+'/*'):
			if filename != description_out and not 'KEGG' in filename:
				with open(filename, 'rb') as readfile:
					shutil.copyfileobj(readfile, outfile)

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
	parser.add_argument('--strat_out', default=False, action='store_true', help='Output table stratified by sequences as well. By ''default this will be in \"contributional\" format ''(i.e. long-format) unless the \"--wide_table\" ''option is set. The startified outfile is named ''\"metagenome_contrib.tsv.gz\" when in long-format.')
	parser.add_argument('--skip_descriptions', default=False, action='store_true', help='Skipping add_descriptions.py step that add a description column to the function abundance table. (default: False')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input_biom', required=True, type=str, help='Input table of sequence abundances (BIOM file used in FPStep1).')
	group_input.add_argument('-f', '--function', required=True, type=str, help='Table of predicted gene family copy numbers ''(output of hsp.py).')
	group_input.add_argument('-m', '--marker', required=True, type=str, help='Table of predicted marker gene copy numbers ''(output of hsp.py - typically for 16S).')
	group_input.add_argument('--max_nsti', type=float, default=2.0, help='Sequences with NSTI values above this value will ' 'be excluded (default: %(default)d).')
	group_input.add_argument('--min_reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for ''each input ASV. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	group_input.add_argument('--min_samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an ASV needs to be ''identfied within. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-o', '--out_dir', metavar='PATH', type=str, default='metagenome_out', help='Output directory for metagenome predictions. ''(default: %(default)s).')
	group_output.add_argument('-t', '--output-tsv', default='FPStep3_abundance.tsv', help='This output file will contain the abundance and metadata (format: TSV). [Default: %(default)s]' )
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.marker)[0])
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		#temp tsv file necessary for metagenome_pipeline.py
		tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
		Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )

		tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )	
		MetagenomePipeline(tmp_biom_to_tsv, args.marker, args.function, args.out_dir, args.max_nsti, args.min_reads, args.min_samples, args.strat_out, tmp_metag_pipeline).submit( args.log_file )

		if not args.skip_descriptions:
			tmp_description_file = tmp_files.add('descriptions_file.tsv.gz')
			formate_description_file(DESCRIPTION_DIR, tmp_description_file )

			pred_file = args.out_dir + "/pred_metagenome_unstrat.tsv.gz"
			AddDescriptions(pred_file,  tmp_description_file, pred_file).submit( args.log_file)
	finally:
		if not args.debug:
			tmp_files.deleteAll()