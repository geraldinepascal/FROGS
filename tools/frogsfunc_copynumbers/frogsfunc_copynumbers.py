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
__version__ = '4.1.0'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'dev'

import os
import re
import sys
import json
import gzip
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
	@summary: Predict number of marker copies (16S, 18S or ITS) for each cluster sequence (i.e OTU).
	"""
	def __init__(self, observed_marker_table, in_tree, hsp_method, output, result_file, log):
		"""
		@param observed_marker_table: [str] Path to marker table file if marker studied is not 16S.
		@param in_tree: [str] Path to resulting tree file with inserted clusters sequences from frogsfunc_placeseqs.
		@param hsp_method: [str] HSP method to use.
		@param output: [str] PICRUSt2 marker output file.
		@param result_file: [str] frogsfunc_copynumbers formatted output file.
		"""
		if observed_marker_table is None:
			input_marker = " -i 16S"
		else:
			input_marker = " --observed_trait_table " + observed_marker_table		

		Cmd.__init__(self,
				 'hsp.py',
				 'predict gene copy number per sequence.', 
				 input_marker + " -t " + in_tree + " --hsp_method " + hsp_method + " -o " + output + " -n  2> " + log,
				"--version")

		self.output = output
		self.result_file = result_file

	def get_version(self):
		return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

	def parser(self, log_file):
		"""
		@summary: Write first column (Clusters names) and last column (nsti score) into final functions results file.
		"""
		if is_gzip(self.output):
			FH_in = gzip.open(self.output,'rt')
			FH_results = gzip.open(self.result_file,'wt')
		else:
			FH_in = open(self.output,'rt')
			FH_results = open(self.result_file,'wt')
	
		for line in FH_in:
			cluster = line.strip().split('\t')[0]
			nsti = line.strip().split('\t')[2]
			FH_results.write(cluster + "\t" + nsti + "\n")
		FH_results.close()

class HspFunction(Cmd):
	"""
	@summary: Predict number of genes family for each cluster sequence (i.e OTU).
	"""
	def __init__(self, in_trait, observed_trait_table, in_tree, hsp_method, output, result_file, log):
		"""
		@param in_trait: [str] Database ID if marker studied is 16S.
		@param observed_trait_table: [str] Path to database trait table if marker studied is not 16S.
		@param in_tree: [str] Path to resulting tree file with insert clusters sequences from frogsfunc_placeseqs.
		@param hsp_method: [str] HSP method to use.
		@param output: [str] PICRUSt2 marker output file.
		@param result_file: [str] frogsfunc_copynumbers formatted output file.
		"""
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
		return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

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

def write_summary(biom_file, output_marker, summary_file ):
	"""
	@summary: Writes the informations to generate graph of the number of OTUs and sequences removed according NSTI score.
	@param biom_file: [str] Biom file (output of frogsfunc_placeseqs).
	@param output_marker: [str] HspMarker step output file, used to find NSTI score per Cluster.
	@param depth_nsti_file: [str] Writes log of nb OTUs and nb sequences kept according to NSTI score.
	"""
	depth_nsti = open(output_marker).readlines()
	max_nsti = 0
	for li in depth_nsti[1:]:
		li = li.strip().split('\t')
		if float(li[2]) > max_nsti:
			max_nsti = float(li[2])
	max_nsti = math.ceil( max_nsti * 50 + 1)

	biom=BiomIO.from_json(biom_file)
	# FH_log = Logger( depth_nsti_file )
	# FH_log.write("#nsti\tnb_clust_kept\tnb_abundances_kept\n")
	step_nsti = [i/50 for i in range(0,max_nsti)] 
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
		# FH_log.write("\t".join([str(nsti), str(len(clusters['Nb'])), str(clusters['Abundances']) ])+"\n")
		
	nstis = sorted(nstis)
	clusters_size = sorted(clusters_size)
	abundances_size = sorted(abundances_size)
	total_abundances = abundances_size[-1]

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "frogsfunc_copynumbers_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )
	for line in FH_summary_tpl:
		if "###CLUSTERS_SIZES###" in line:
			line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
		elif "###ABUNDANCES_SIZES###" in line:
			line = line.replace( "###ABUNDANCES_SIZES###", json.dumps( abundances_size) )
		elif "###NSTI_THRESH###" in line:
			line = line.replace( "###NSTI_THRESH###", json.dumps(nstis) )
		FH_summary_out.write( line )
	# FH_log.close()
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
	group_input.add_argument('-b', '--input-biom', required=True, help='frogsfunc_placeseqs output Biom file (frogsfunc_placeseqs.biom).')
	group_input.add_argument('-t', '--input-tree', required=True, type=str, help='frogsfunc_placeseqs output tree in newick format containing both studied sequences (i.e. ASVs or OTUs) and reference sequences.')
	group_input.add_argument('--hsp-method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use. mp: predict discrete traits using max parsimony. emp_prob: predict discrete traits based on empirical state probabilities across tips. subtree_average: predict continuous traits using subtree averaging. pic: predict continuous traits with phylogentic independent contrast. scp: reconstruct continuous traits using squared-change parsimony (default: %(default)s).')

	group_input.add_argument('-m', '--marker-type', required=True, choices=['16S','ITS','18S'], help='Marker gene to be analyzed.')
	group_input_16S = parser.add_argument_group( '16S ' )
	group_input_16S.add_argument('-i', '--input-functions', default=["EC"], nargs='+', choices=['EC', 'KO', 'COG', 'PFAM', 'TIGRFAM','PHENO'], help="Specifies which function databases should be used (%(default)s). EC is used by default because necessary for frogsfunc_pathways. At least EC or KO is required. To run the command with several functions, separate the functions with spaces (ex: -i EC PFAM).")

	group_input_other = parser.add_argument_group( 'ITS and 18S ' )
	group_input_other.add_argument('--input-function-table',help="The path to input functions table describing directly observed functions, in tab-delimited format.(ex $PICRUSt2_PATH/default_files/fungi/ec_ITS_counts.txt.gz). Required.")
	group_input_other.add_argument('--input-marker-table',help="The input marker table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. (ex $PICRUSt2_PATH/default_files/fungi/ITS_counts.txt.gz). Required.")

	# Output
	group_output = parser.add_argument_group( 'Outputs' )
	group_output.add_argument('-o', '--output-marker', default="frogsfunc_copynumbers_marker.tsv", type=str, help='Output table of predicted marker gene copy numbers per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.[Default: %(default)s]')
	group_output.add_argument('-f', '--output-function', default="frogsfunc_copynumbers_predicted_functions.tsv", type=str, help='Output table with predicted function abundances per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.[Default: %(default)s]')
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-s', '--summary', default='frogsfunc_copynumbers_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)

	# Check for 16S input
	if args.marker_type == "16S" and (not 'EC' in args.input_functions and not 'KO' in args.input_functions):
			parser.error("\n\n#ERROR : --input-functions : 'EC' and/or 'KO' must be at least indicated (others functions are optionnal)")
	# Check for ITS or 18S input
	if args.marker_type in ["ITS", "18S"]:
		if (args.input_function_table is not None and args.input_marker_table is None) or (args.input_function_table is None and args.input_marker_table is not None):
			parser.error("\n\n#ERROR : --input-function-table and --input-marker-table both required when studied marker is not 16S!\n\n")
		elif args.input_function_table is not None and args.input_marker_table is not None:
			args.input_functions = None

	tmp_files=TmpFiles(os.path.split(args.output_marker)[0])

	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		tmp_hsp_marker = tmp_files.add( 'tmp_hsp_marker.log' )
		HspMarker(args.input_marker_table, args.input_tree, args.hsp_method, args.output_marker, args.output_function, tmp_hsp_marker).submit(args.log_file)

		# tmp_depth_nsti = tmp_files.add( 'depth_nsti.txt' )
		write_summary(args.input_biom, args.output_marker, args.summary)

		tmp_hsp_function = tmp_files.add( 'tmp_hsp_function.log' )
		if args.input_functions is not None:

			suffix_name = "_copynumbers_predicted.tsv"
			for trait in args.input_functions:
				cur_output_function = trait + suffix_name
				tmp_output_function = tmp_files.add( cur_output_function )
				Logger.static_write(args.log_file, '\n\nRunning ' + trait + ' functions prediction.\n')
				HspFunction(trait, args.input_function_table, args.input_tree, args.hsp_method, tmp_output_function, args.output_function, tmp_hsp_function).submit(args.log_file)
		else:
			cur_output_function = "copynumbers_predicted.tsv"
			tmp_output_function = tmp_files.add( cur_output_function )
			HspFunction(args.input_functions, args.input_function_table, args.input_tree, args.hsp_method, tmp_output_function, args.output_function, tmp_hsp_function).submit(args.log_file)
	
	finally:
		if not args.debug:
			tmp_files.deleteAll()
