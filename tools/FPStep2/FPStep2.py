#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard & Vincent Darbot & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import re
import os
import sys
import json
import gzip
import argparse
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
	@summary: hsp.py : predict number of marker copies (16S, 18S or ITS) for each OTU.

	"""
	def __init__(self, observed_marker_table, in_tree, hsp_method, output, result_file, log):
		if observed_marker_table is None:
			input_marker = " -i 16S"
		else:
			input_marker = " --observed_trait_table " + observed_marker_table		

		Cmd.__init__(self,
				 'hsp.py',
				 'predict gene copy numer per sequence.', 
				 input_marker + " -t " + in_tree + " --hsp_method " + hsp_method + " -o " + output + " -n  2> " + log,
				"--version")

		self.output = output
		self.result_file = result_file

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

	def parser(self, log_file):
		"""
		@summary: Write first column (Clusters names) and last column (nsti score) into final functions results file.
		"""
		if is_gzip(self.output):
			FH_in = gzip.open(self.output,'rt').readlines()
			FH_results = gzip.open(self.result_file,'wt')
		else:
			FH_in = open(self.output,'rt').readlines()
			FH_results = open(self.result_file,'wt')
	
		for line in FH_in:
			cluster = line.strip().split('\t')[0]
			nsti = line.strip().split('\t')[2]
			FH_results.write(cluster + "\t" + nsti + "\n")
		FH_results.close()

class HspFunction(Cmd):
	"""
	@summary: hsp.py predict number of genes family for each OTU.
	"""
	def __init__(self, in_trait, observed_trait_table, in_tree, hsp_method, output, result_file, log):
		if observed_trait_table is None:
			input_function = " --in_trait " + in_trait
		else:
			input_function = " --observed_trait_table " + observed_trait_table

		Cmd.__init__(self,
				 'hsp.py',
				 'predict function abundance per sequence.', 
				  input_function + " -t " + in_tree + " --hsp_method " + hsp_method +  " -o " + output + " -n 2>> " + log,
				"--version")

		self.output = output
		self.result_file = result_file

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

	def parser(self, log_file):
		"""
		@summary: Concatane function tables of predicted abundances into one global.
		"""
		if is_gzip(self.output):
			FH_in = gzip.open(self.output,'rt').readlines()
			tmp = gzip.open(self.result_file+'.tmp', 'wt')
			FH_results = gzip.open(self.result_file,'rt').readlines()
		else:
			FH_in = open(self.output,'rt').readlines()
			tmp = open(self.result_file+'.tmp', 'wt')
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

def is_gzip( file ):
	"""
	@return: [bool] True if the file is gziped.
	@param file: [str] Path to processed file.
	"""
	is_gzip = None
	FH_input = gzip.open( file )
	try:
		FH_input.readline()
		is_gzip = True
	except:
		is_gzip = False
	finally:
		FH_input.close()
	return is_gzip

def check_functions( functions ):
	"""
	@summary: check if --function parameter is valid.
	"""
	VALID_FUNCTIONS = ['EC','COG','KO','PFAM','TIGRFAM','PHENO']
	# if the user add mulitple functions prediction
	if "," in functions:
		functions = functions.split(',')
	else:
		functions = [functions]
	for function in functions:
		if function not in VALID_FUNCTIONS:
			raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : With '--function' parameter: " + function + " not a valid function.\n\n" ))
	return functions

def write_summary(biom_file, output_marker, depth_nsti_file, summary_file ):
	"""
	"""
	depth_nsti = open(output_marker).readlines()
	biom=BiomIO.from_json(biom_file)
	FH_log = Logger( depth_nsti_file )
	FH_log.write("#nsti\tnb_clust_kept\tnb_abundances_kept\n")
	step_nsti = [i/50 for i in range(0,101)] 
	cluster_kept = dict()
	for cur_nsti in step_nsti:
		cluster_kept[cur_nsti] = { 'Nb' : [], 'Abundances' : 0 }
		for li in depth_nsti[1:]:
			li = li.strip().split('\t')
			if float(li[2]) <= cur_nsti:
				cluster_kept[cur_nsti]['Nb'].append(li[0])
				cluster_kept[cur_nsti]['Abundances']+=biom.get_observation_count(li[0])
	clusters_size = list()
	abundances_size = list()
	nstis = list()
	for nsti,clusters in cluster_kept.items():
		clusters_size.append(len(clusters['Nb']))
		abundances_size.append(clusters['Abundances'])
		nstis.append(float(nsti))
		FH_log.write("\t".join([str(nsti), str(len(clusters['Nb'])), str(clusters['Abundances']) ])+"\n")
	nstis = sorted(nstis)
	clusters_size = sorted(nstis)
	abundances_size = sorted(abundances_size)
	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep2_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )
	for line in FH_summary_tpl:
		if "###CLUSTERS_SIZES###" in line:
			line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
		elif "###ABUNDANCES_SIZES###" in line:
			line = line.replace( "###ABUNDANCES_SIZES###", json.dumps( abundances_size) )
		elif "###NSTI_THRESH###" in line:
			line = line.replace( "###NSTI_THRESH###", json.dumps(nstis) )
		FH_summary_out.write( line )
	FH_log.close()
	FH_summary_out.close()
	FH_summary_tpl.close()

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser( description='predict gene family for OTU' )
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input_biom', required=True, help='Biom file.')
	group_input.add_argument('-i', '--in_trait', default="EC",help="For 16S marker input: Specifies which default trait table should be used ('EC', 'KO', 'COG', PFAM', 'TIGRFAM' or 'PHENO'). EC is used by default because necessary for FPStep4. To run the command with several functions, separate the functions with commas (ex: KO,PFAM). (for ITS or 18S : only EC available)")
	group_input.add_argument('--observed_marker_table',help="The input marker table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you don't work on 16S marker. (ex $PICRUST2_PATH/default_files/fungi/ITS_counts.txt.gz). This input is required when the --observed_trait_table option is set. ")
	group_input.add_argument('--observed_trait_table',help="The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you don't work on 16S marker. (ex $PICRUST2_PATH/default_files/fungi/ec_ITS_counts.txt.gz). This input is required when the --observed_marker_table option is set. ")
	group_input.add_argument('-t', '--tree', required=True, type=str, help='FPStep1 output tree in newick format containing both study sequences (i.e. ASVs or OTUs) and reference sequences.')
	group_input.add_argument('-s', '--hsp_method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use.' +'"mp": predict discrete traits using max parsimony. ''"emp_prob": predict discrete traits based on empirical ''state probabilities across tips. "subtree_average": ''predict continuous traits using subtree averaging. ' '"pic": predict continuous traits with phylogentic ' 'independent contrast. "scp": reconstruct continuous ''traits using squared-change parsimony (default: ''%(default)s).')
	# Output
	group_output = parser.add_argument_group( 'Outputs' )
	group_output.add_argument('-m', '--output_marker', default="FPStep2_marker_copy_number.tsv", type=str, help='Output table of predicted marker gene copy numbers per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')
	group_output.add_argument('-o', '--output_function', default="FPStep2_predicted_functions.tsv", type=str, help='Output table with predicted abundances per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('--html', default='FPStep2_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)

	if args.in_trait is not None:
		if not 'EC' in args.in_trait and not 'KO' in args.in_trait:
			parser.error("\n\n#ERROR : --in_trait : 'EC' and/or 'KO' must be at least indicated (others functions are optionnal)")
	if (args.observed_trait_table is not None and args.observed_marker_table is None) or (args.observed_trait_table is None and args.observed_marker_table is not None):
		parser.error("\n\n#ERROR : --observed_trait_table and --observed_marker_table both required when studied marker is not 16S!\n\n")
	elif args.observed_trait_table is not None and args.observed_marker_table is not None:
		args.in_trait = None




	# default output marker file name
	if args.output_marker is None:
		args.output_marker = "FPStep2_marker_nsti_predicted.tsv"

	tmp_files=TmpFiles(os.path.split(args.output_marker)[0])

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		tmp_hsp_marker = tmp_files.add( 'tmp_hsp_marker.log' )
		HspMarker(args.observed_marker_table, args.tree, args.hsp_method, args.output_marker, args.output_function, tmp_hsp_marker).submit(args.log_file)

		tmp_depth_nsti = tmp_files.add( 'depth_nsti.txt' )
		write_summary(args.input_biom, args.output_marker, tmp_depth_nsti, args.html)

		tmp_hsp_function = tmp_files.add( 'tmp_hsp_function.log' )
		if args.in_trait is not None:
			in_traits = check_functions(args.in_trait)

			suffix_name = "_predicted.tsv"
			for trait in in_traits:
				cur_output_function = trait + suffix_name
				Logger.static_write(args.log_file, '\n\nRunning ' + trait + ' functions prediction.\n')
				HspFunction(trait, args.observed_trait_table, args.tree, args.hsp_method, cur_output_function, args.output_function, tmp_hsp_function).submit(args.log_file)
		else:
			cur_output_function = "function_predicted.tsv"
			HspFunction(args.in_trait, args.observed_trait_table, args.tree, args.hsp_method, cur_output_function, args.output_function, tmp_hsp_function).submit(args.log_file)
	
	finally:
		if not args.debug:
			tmp_files.deleteAll()