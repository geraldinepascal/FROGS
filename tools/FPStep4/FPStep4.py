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
from numpy import median
import gzip
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

#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

# import picrust2
from picrust2.pathway_pipeline import pathway_pipeline
from picrust2.util import (make_output_dir, check_files_exist,
						   TemporaryDirectory)
from picrust2.default import default_regroup_map, default_pathway_map
from os import path
##################################################################################################################################################
#
# COMMAND LINES #  pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \ -o pathways_out \ --intermediate minpath_working \ -p 1
				
##################################################################################################################################################
class pathway_pipeline(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: predict gene abudance
	@see: https://github.com/picrust/picrust2/wiki
	"""
	
	def __init__(self, pathwayMet, input_file, output, stdout):
	
		Cmd.__init__(self,
				 'pathway_pipeline.py ',
				 'predict abundance pathway', 
				  "" + pathwayMet  +" -i "+ str(input_file) +" -o "+ str(output) +' 2> ' + stdout,
				"--version") 
	  
		
	# def get_version(self):
	#	 """
	#	 @summary: Returns the program version number.
	#	 @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
	#	 """

	#	 return Cmd.get_version(self, 'stdout').split()[1].strip() 

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
"""
			derive of picrust2 : Start
"""

def pathway_pipeline(inputfile,
					 mapfile,
					 out_dir,
					 proc=1,
					 run_minpath=True,
					 coverage=False,
					 no_regroup=False,
					 regroup_mapfile=None,
					 gap_fill_on=True,
					 per_sequence_contrib=False,
					 per_sequence_abun=None,
					 per_sequence_function=None,
					 wide_table=False,
					 verbose=False):
	'''Pipeline containing full pipeline for reading input files, making
	calls to functions to run MinPath and calculate pathway abundances and
	coverages. Will return: (1) unstratified pathway abundances, (2)
	unstratified pathway coverages, (3) stratified pathway abundances, (4)
	stratified pathway coverages, (5) pathway abundance predictions per
	sequence, (6) pathway coverage predictions per sequence, (7) unstratified
	pathway abundances based on per-sequence predictions. An object of class
	None will be returned for any non-applicable value.'''

	# If no regrouping flag set then set input regrouping mapfile to be None.
	if no_regroup:
		regroup_mapfile = None

	# Read in table of gene family abundances and determine if in unstratified,
	# stratified, or contribution format.
	in_metagenome, in_format = read_metagenome_input(inputfile)

	# Basic checks if --per_sequence_contrib set.
	if per_sequence_contrib:

		# Throw error if --per_sequence_contrib set, but --per_sequence_abun
		# and/or --per_sequence_function not set.
		if not per_sequence_abun or not per_sequence_function:
			sys.exit("Error: \"--per_sequence_contrib\" option set, but at "
					 "least one of \"per_sequence_abun\" or "
					 "\"--per_sequence_function\" were not set. These input "
					 "arguments need to be specified when "
					 "\"--per_sequence_contrib\" is used.")

		check_files_exist([per_sequence_abun, per_sequence_function])

	# Throw error file-format and wide table setting not compatible.
	if in_format == "strat" and not wide_table:
		sys.exit("Error: stratified table input (deprecated format), but "
				 "\"--wide_table\" option not set. You should input either an "
				 "unstratified or contributional table if you do not require "
				 "a wide-format table.")

	if in_format == "contrib" and wide_table and not per_sequence_contrib:
		sys.exit("Error: contributional table input, but \"--wide_table\" "
				 "option set. This option specifies that deprecated "
				 "wide-format stratified tables should be output, which "
				 "is only allowed when a wide-format stratified table is "
				 "input or the --per_sequence_contrib option is set.")

	# Remove 'description' column if it exists.
	if "description" in in_metagenome.columns:
		in_metagenome.drop("description", axis=1, inplace=True)

	# Get list of sample ids.
	if in_format == "contrib":
		samples = in_metagenome['sample'].unique()
	else:
		samples = [col for col in in_metagenome.columns
				   if col not in ["function", "sequence"]]

	# Initialize reactions to be empty unless regroup mapfile given.
	reactions = []

	# Regroup functions in input table to different ids if regroup mapfile is
	# provided.
	if regroup_mapfile:
		reactions = read_reaction_names(regroup_mapfile)

		in_metagenome = regroup_func_ids(in_metagenome, in_format,
										 regroup_mapfile, proc)
		regrouped_outfile = path.join(out_dir, "regrouped_infile.tsv")
		in_metagenome.to_csv(path_or_buf=regrouped_outfile, sep="\t",
							 index=False)

	# Read in pathway structures.
	pathways_in = PathwaysDatabase(database=mapfile, reaction_names=reactions)

	# Write out mapfile with all structure removed.
	if run_minpath:
		minpath_mapfile = path.join(out_dir, "parsed_mapfile.tsv")
		with open(minpath_mapfile, "w") as out_map:
			out_map.write(pathways_in.get_database())
	else:
		minpath_mapfile = None

	# Subset input table of reactions to only those found in pathway database.
	in_metagenome = in_metagenome[in_metagenome.function.isin(pathways_in.reaction_list())]

	# Initialize output objects to be None (expect for unstratified abundance).
	path_cov_unstrat = None
	path_cov_strat = None
	path_abun_strat = None
	path_cov_by_seq = None
	path_abun_by_seq = None
	path_abun_unstrat_by_seq = None

	minpath_out_dir = path.join(out_dir, "minpath_running")
	make_output_dir(minpath_out_dir)

	if in_format == "contrib":
		# Get unstratified and stratified pathway levels.
		# Note that stratified tables will only be returned by this step (and
		# the "strat" option below) if per_sequence_contrib=False (extra step
		# required below).
		path_out_raw = Parallel(n_jobs=proc)(delayed(contrib_format_pathway_levels)(
													 sample_id,
													 in_metagenome.loc[in_metagenome['sample'] == sample_id],
													 minpath_mapfile,
													 minpath_out_dir,
													 pathways_in,
													 run_minpath,
													 coverage,
													 gap_fill_on,
													 per_sequence_contrib,
													 verbose)
													 for sample_id in samples)

	elif in_format == "strat":

		path_out_raw = Parallel(n_jobs=proc)(delayed(basic_strat_pathway_levels)(
													 sample_id,
													 in_metagenome[["function", "sequence", sample_id]],
													 minpath_mapfile,
													 minpath_out_dir,
													 pathways_in,
													 run_minpath,
													 coverage,
													 gap_fill_on,
													 per_sequence_contrib,
													 verbose)
													 for sample_id in samples)

	# Otherwise the data is in unstratified format, which is more straight-
	# forward to process.
	else:
		path_out_raw = Parallel(n_jobs=proc)(delayed(
											   unstrat_pathway_levels)(
												   sample_id,
												   in_metagenome[["function", sample_id]],
												   minpath_mapfile,
												   minpath_out_dir,
												   pathways_in,
												   run_minpath,
												   coverage,
												   gap_fill_on,
												   verbose)
											   for sample_id in samples)

	# Prep output unstratified DataFrames.
	path_raw_abun_unstrat = []
	path_raw_cov_unstrat = []

	for sample_output in path_out_raw:
		path_raw_abun_unstrat += [sample_output[0]]
		path_raw_cov_unstrat += [sample_output[1]]

	path_abun_unstrat = prep_pathway_df_out(path_raw_abun_unstrat)
	path_abun_unstrat.columns = samples
	path_abun_unstrat.sort_index(axis=0, inplace=True)

	if coverage:
		path_cov_unstrat = prep_pathway_df_out(path_raw_cov_unstrat,
											   num_digits=10)
		path_cov_unstrat.columns = samples
		path_cov_unstrat.sort_index(axis=0, inplace=True)
	else:
		path_cov_unstrat = None

	# If --per_sequence_contrib not set then prep output stratified
	# table the same as the unstratified tables.
	if not per_sequence_contrib and in_format != "unstrat":
		path_raw_abun_strat = []

		for sample_output in path_out_raw:
			path_raw_abun_strat += [sample_output[2]]

		if in_format == "strat":
			path_abun_strat = prep_pathway_df_out(path_raw_abun_strat,
												  strat_index=True)
			path_abun_strat.columns = ["pathway", "sequence"] + samples
			path_abun_strat.sort_values(['pathway', 'sequence'], inplace=True)

		elif in_format == "contrib":
			path_abun_strat = pd.concat(path_raw_abun_strat)
			path_abun_strat.sort_values(['sample', 'function', 'taxon'],
										inplace=True)

	# Calculate pathway levels for each individual sequence (in parallel)
	# and then multiply this table by the abundance of each sequence
	# within each sample (using same approach as in metagenome pipeline).
	if per_sequence_contrib:

		per_seq_out_dir = path.join(out_dir, "minpath_running_per_seq")
		make_output_dir(per_seq_out_dir)

		path_abun_strat, \
		path_cov_strat, \
		path_abun_by_seq, \
		path_cov_by_seq = per_sequence_contrib_levels(sequence_abun=per_sequence_abun,
													  sequence_func=per_sequence_function,
													  minpath_map=minpath_mapfile,
													  per_seq_out_dir=per_seq_out_dir,
													  pathway_db=pathways_in,
													  run_minpath=run_minpath,
													  calc_coverage=coverage,
													  gap_fill_on=gap_fill_on,
													  nproc=proc,
													  regroup_map=regroup_mapfile,
													  wide_table=wide_table,
													  print_opt=verbose)

		if wide_table:
			path_abun_unstrat_by_seq = strat_to_unstrat_counts(strat_df=path_abun_strat,
															   func_col="pathway")
		else:
			path_abun_unstrat_by_seq = contrib_to_unstrat(contrib_table=path_abun_strat,
														  sample_order=list(path_abun_unstrat.columns.values))

	return(path_abun_unstrat, path_cov_unstrat, path_abun_strat,
		   path_cov_strat, path_abun_by_seq, path_cov_by_seq,
		   path_abun_unstrat_by_seq)

def make_output_dir(dirpath, strict=False):
	"""Make an output directory if it doesn't exist

	Returns the path to the directory
	dirpath -- a string describing the path to the directory
	strict -- if True, raise an exception if dir already
	exists
	"""
	dirpath = abspath(dirpath)

	# Check if directory already exists
	if isdir(dirpath):
		if strict:
			err_str = "Directory '%s' already exists" % dirpath
			raise IOError(err_str)

		return dirpath
	try:
		makedirs(dirpath)
	except IOError as e:
		err_str = "Could not create directory '%s'. Are permissions set " +\
				  "correctly? Got error: '%s'" %e
		raise IOError(err_str)

	return dirpath

def check_files_exist(filepaths):
	'''Takes in a list of filepaths and checks whether they exist. Will
	throw error describing which files do not exist if applicable.'''

	num_nonexist = 0

	missing_files = []

	for filepath in filepaths:

		if not exists(filepath):
			missing_files += [filepath]
			num_nonexist += 1

	if num_nonexist == 0:
		pass
	elif num_nonexist == 1:
		raise ValueError("This input file was not found: " + missing_files[0])
	elif num_nonexist > 1:
		raise ValueError("These input files were not found: " +
						 ", ".join(missing_files))

class TemporaryDirectory(object):
	'''Create and return a temporary directory.  This has the same
	behavior as mkdtemp but can be used as a context manager.  For
	example:
		with TemporaryDirectory() as tmpdir:
			...
	Upon exiting the context, the directory and everything contained
	in it are removed.

	NOTE: This function was taken and modified from the tempfile package to
	first change permissions on folder to be deleted.'''

	def __init__(self, suffix=None, prefix=None, dir=None):
		self.name = tempfile.mkdtemp(suffix, prefix, dir)
		self._finalizer = _weakref.finalize(
			self, self._cleanup, self.name,
			warn_message="Implicitly cleaning up {!r}".format(self))

	@classmethod
	def _cleanup(cls, name, warn_message):
		_shutil.rmtree(name)
		_warnings.warn(warn_message, ResourceWarning)

	def __repr__(self):
		return "<{} {!r}>".format(self.__class__.__name__, self.name)

	def __enter__(self):
		return self.name

	def __exit__(self, exc, value, tb):
		self.cleanup()

	def cleanup(self):
		if self._finalizer.detach():

			# Line added by Gavin Douglas to change permissions to 777 before
			# deleting:
			call(["chmod", "-R", "777", self.name])

			_shutil.rmtree(self.name)

"""
			derive of picrust2 : End
	
"""

def dezip(file1, out_file1):

	#file_zip = gzip.GzipFile("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv.gz", 'rb')
	file_zip = gzip.GzipFile(file1, 'rb')
	s = file_zip.read()
	file_zip.close()

	#file_dezip = open ("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv", 'wb')
	file_dezip = open (out_file1, 'wb')
	file_dezip.write(s)
	file_dezip.close()

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":

	# Manage parameters
	parser = argparse.ArgumentParser( description='Phylogenetic tree reconstruction' )
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )

	group_input.add_argument('-i', '--input_file', metavar='IN_TABLE', required=True, type=str, help='Input TSV table of gene family abundances (either ''the unstratified or stratified output of ' 'metagenome_pipeline.py).')
	

	#Outputs
	group_output = parser.add_argument_group( 'Outputs')

	#group_output.add_argument('-o', '--output', metavar='DIRECTORY', required=True, type=str, help='Output folder for pathway abundance output.')

	group_output.add_argument('--path_abund', metavar='PATH', default='path_abun_unstrat.tsv.gz', help='Output file for metagenome predictions abundance. ''(default: %(default)s).')
	
	group_output.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)

	group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

	args = parser.parse_args()

	stderr = "FPStep4.stderr"

	pathwayMet = ""
	# Process 
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		 ########### Création de chemin ############# 
		output = os.path.abspath(os.path.dirname(args.path_abund)) + "/" +str(time.time()) + "_" + str(os.getpid())
		###### Lancer commande #####
		pathway_pipeline(pathwayMet, args.input_file, output, stderr).submit(args.log_file)
		####### Mettre le fichier dans le nouveau dossier (de l'argument path_abund)
		os.system("mv " + output +"/"+ "path_abun_unstrat.tsv.gz " + args.path_abund)
		# Dezip la sortie zipper
		dezip(args.path_abund, "path_abun_unstrat.tsv")
	finally:
		print("Partie Finale ")


