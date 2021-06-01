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
import sys
import argparse
import json
import re
import gzip
from collections import OrderedDict

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPAT
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
# Default table PATH
ITS_PATH = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/ITS_counts.txt.gz")
ITS_EC_PATH = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/ec_ITS_counts.txt.gz")
_18S_PATH = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/18S_counts.txt.gz")
_18S_EC_PATH = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/ec_18S_counts.txt.gz")
#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class HspMarker(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: hsp.py program  predict number copie of gene family for each OTU.
	@summary: use 16S, EC and/or KO 
	@see: https://github.com/picrust/picrust2/wiki
	"""
	def __init__(self, category, in_tree, output, log):
		if category == "16S":
			input_marker = " -i 16S"
		elif category == "ITS":
			input_marker = " --observed_trait_table " + ITS_PATH
		elif category == "18S":
			input_marker = " --observed_trait_table " + _18S_PATH			

		Cmd.__init__(self,
				 'hsp.py',
				 'predict gene copy numer per sequence.', 
				 input_marker + " -t " + in_tree + " -o " + output + " -n  2> " + log,
				"--version")

	def get_version(self):
		"""
		@summary: Returns the program version number.
		@return: [str] Version number if this is possible, otherwise this method return 'unknown'.
		"""
		return Cmd.get_version(self, 'stdout').split()[1].strip()

class HspFunction(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: hsp.py program  predict number copie of gene family for each OTU.
	@summary: use 16S, EC and/or KO 
	@see: https://github.com/picrust/picrust2/wiki
	"""
	def __init__(self, category, function, in_tree, output, result_file, log):
		if category == "16S":
			input_function = " -i " + function
		elif category == "ITS":
			input_function = " --observed_trait_table " + ITS_EC_PATH
		elif category == "18S":
			input_function = " --observed_trait_table " + _18S_EC_PATH
		Cmd.__init__(self,
				 'hsp.py',
				 'predict function abundance per sequence.', 
				  input_function + " -t " + in_tree + " -o " + output + " -n 2>> " + log,
				"--version")

		self.output = output
		self.result_file = result_file

	def get_version(self):
		"""
		@summary: Returns the program version number.
		@return: [str] Version number if this is possible, otherwise this method return 'unknown'.
		"""
		return Cmd.get_version(self, 'stdout').split()[1].strip()

	def parser(self, log_file):
		"""
		@summary: Concatane function tables of predicted abundances into one global.
		"""
		FH_in = gzip.open(self.output,'rt').readlines()
		# if it's the first function
		if not os.path.exists(self.result_file):
			FH_results = open(self.result_file,'w')
			for line in FH_in:
				FH_results.write(line)
			FH_results.close()
		# we add columns of others tables without the first column (Cluster names) and last (nsti score).
		else:
			tmp = open(self.result_file+'.tmp', 'w')
			FH_results = open(self.result_file,'rt').readlines()
			for cur_line in range(len(FH_in)):
				line = FH_in[cur_line].split('\t')[1:-1]
				result = FH_results[cur_line].split('\t')
				tmp.write("\t".join(result[0:-1])+"\t"+"\t".join(line)+"\t"+result[-1])
			tmp.close()
			os.rename(self.result_file+'.tmp', self.result_file)

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser( description='predict gene family for OTU' )
	parser.add_argument('--verbose', default=False, action='store_true',help='If specified, print out wrapped commands and other ''details to screen.')
	parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-c', '--category', choices=['16S', '18S', 'ITS'], default='16S', help='Specifies which category 16S, 18S or ITS')
	group_input.add_argument('-f', '--function', default=None,help="Specifies which default trait table should be used ('EC', 'COG', 'KO', 'PFAM', 'TIGRFAM' or 'PHENO'). To run the command with several functions, separate the functions with commas (ex: COG,KO,PFAM) (for ITS or 18S : only EC available")
	group_input.add_argument('-t', '--tree', required=True, type=str, help='FPStep1 output tree in newick format containing both study sequences (i.e. ASVs or OTUs) and reference sequences.')
	group_input.add_argument('-s', '--hsp_method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use.' +'"mp": predict discrete traits using max parsimony. ''"emp_prob": predict discrete traits based on empirical ''state probabilities across tips. "subtree_average": ''predict continuous traits using subtree averaging. ' '"pic": predict continuous traits with phylogentic ' 'independent contrast. "scp": reconstruct continuous ''traits using squared-change parsimony (default: ''%(default)s).')
	# Output
	group_output = parser.add_argument_group( 'Outputs' )
	group_output.add_argument('-m', '--output_marker', default=None, type=str, help='Output table of predicted marker gene copy numbers per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')
	group_output.add_argument('-o', '--output_function', default=None, type=str, help='Output table with predicted abundances per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	args = parser.parse_args()
	prevent_shell_injections(args)

	# manage category and function input parameters
	if args.category == "16S" and args.function is None:
		args.function = "EC"
	if args.category == "ITS" and args.function is None:
		args.function = "EC"
	if args.category == "ITS" and args.function != "EC":
		Logger.static_write(args.log_file, '\n\n#WARNING : --function parameter: only EC available with --category set to ITS.\n')
		Logger.static_write(args.log_file, '\n--function parameter set to EC.\n\n')
		args.function = "EC"
	# default output marker file name
	if args.output_marker is None:
		args.output_marker = args.category + "_nsti_predicted.tsv.gz"

	tmp_files=TmpFiles(os.path.split(args.output_marker)[0])
	functions_files = []
	# if the user add mulitple functions prediction
	if "," in args.function:
		functions = args.function.split(',')
	else:
		functions = [args.function]

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		tmp_hsp_marker = tmp_files.add( 'tmp_hsp_marker.log' )
		HspMarker(args.category, args.tree, args.output_marker, tmp_hsp_marker).submit(args.log_file)

		tmp_hsp_function = tmp_files.add( 'tmp_hsp_function.log' )

		for function in functions:
			#default output names 
			if args.output_function is None:
				suffix = "_predicted.tsv.gz"
				final_predicted = "all_predicted.tsv.gz"
				cur_output_function = function + suffix
			else:
				suffix = args.output_function
				final_predicted = "all_" + args.output_function
				cur_output_function = function + "_" + suffix
			functions_files.append(cur_output_function)
			Logger.static_write(args.log_file, '\n\nRunning ' + function + ' functions prediction.\n')
			HspFunction(args.category, function, args.tree, cur_output_function, final_predicted, tmp_hsp_function).submit(args.log_file)

	finally:
		if not args.debug:
			tmp_files.deleteAll()