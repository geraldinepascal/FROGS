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
import json
import gzip
import argparse
import pandas as pd
from numpy import median
from shutil import rmtree
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
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: predict gene abudance
	@see: https://github.com/picrust/picrust2/wiki
	"""
	def __init__(self, in_biom, marker, function, out_dir, log):
	
		Cmd.__init__(self,
				 'metagenome_pipeline.py ',
				 'Per-sample functional profiles prediction.', 
				 " -i " +  in_biom + " -m " + marker + " -f " + function + " -o " + out_dir + ' 2> ' + log,
				"--version") 
	  
	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

class Biom2tsv(Cmd):
	"""
	@summary: Converts BIOM file to TSV file.
	@note: taxonomyRDP seedID seedSequence blastSubject blastEvalue blastLength blastPercentCoverage blastPercentIdentity blastTaxonomy OTUname SommeCount sample_count
	"""
	def __init__(self, in_biom, out_tsv):

		Cmd.__init__( self,
					  'biom2tsv.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-file " + in_biom + " --output-file " + out_tsv,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip() 


class Biom2multiAffi(Cmd):
	"""
	@summary: Extracts multi-affiliations from a FROGS BIOM file.
	"""
	def __init__( self, out_tsv, in_biom ):
		"""
		@param in_biom: [str] Path to BIOM file.
		@param out_tsv: [str] Path to output TSV file.
		"""
		# Set command
		Cmd.__init__( self,
					  'multiAffiFromBiom.py',
					  'Extracts multi-affiliations data from a FROGS BIOM file.',
					  '--input-file ' + in_biom + ' --output-file ' + out_tsv,
					  '--version' )

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def tsvParse(out_tsv):

	tsv_parse = "out_tsv.tsv"
	# je crée un dataframe
	df = pd.read_csv(out_tsv, header=0, delimiter='\t')
	# Je suprime les mauvaises colonnes
	al = df.drop(df.columns[[0, 1, 3]], axis=1)
	#print(df.drop(df.columns[[0, 1]], axis=1))
	print(al)
	# J écris dans un nouveau fichier
	al.to_csv(tsv_parse, "\t", index=None)

##### Fonction dezip
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
	parser = argparse.ArgumentParser( description='Per-sample functional profiles prediction.' )
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	parser.add_argument('--strat_out', default=False, action='store_true', help='Output table stratified by sequences as well. By ''default this will be in \"contributional\" format ''(i.e. long-format) unless the \"--wide_table\" ''option is set. The startified outfile is named ''\"metagenome_contrib.tsv.gz\" when in long-format.')
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

		tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
		Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )

		tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )	
		MetagenomePipeline(tmp_biom_to_tsv, args.marker, args.function, args.out_dir, tmp_metag_pipeline).submit( args.log_file )

	finally:
		if not args.debug:
			tmp_files.deleteAll()