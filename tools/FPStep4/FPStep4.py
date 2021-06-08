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

##################################################################################################################################################
#
# COMMAND LINES #  pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \ -o pathways_out \ --intermediate minpath_working \ -p 1
				
##################################################################################################################################################
class PathwayPipeline(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: predict gene abudance
	@see: https://github.com/picrust/picrust2/wiki
	"""
	def __init__(self, input_file, out_dir, log):
	
		Cmd.__init__(self,
				 'pathway_pipeline.py ',
				 'predict abundance pathway', 
				  " -i " + input_file + " -o " + out_dir + ' 2> ' + log,
				"--version") 
	  
		
	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

class AddDescriptions(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: predict gene abudance
	@see: https://github.com/picrust/picrust2/wiki
	"""
	def __init__(self, input_file, out_dir, log):
	
		Cmd.__init__(self,
				 'pathway_pipeline.py ',
				 'predict abundance pathway', 
				  " -i " + input_file + " -o " + out_dir + ' 2> ' + log,
				"--version") 
	  
		
	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

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
	parser = argparse.ArgumentParser( description='Infer the presence and abundances of pathways based on gene family abundances in a sample.' )
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-i', '--input_file', required=True, type=str, help='Input TSV table of gene family abundances (either ''the unstratified or stratified output of ' 'metagenome_pipeline.py).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-o', '--out_dir', default='pathways_out', help='Output folder for pathway abundance output.')
	group_output.add_argument('--path_abund', default='path_abun_unstrat.tsv.gz', help='Output file for metagenome predictions abundance. (default: %(default)s).')
	group_output.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmp_files=TmpFiles(os.path.split(args.input_file)[0])
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		 ########### Création de chemin ############# 
		output = os.path.abspath(os.path.dirname(args.path_abund)) + "/" +str(time.time()) + "_" + str(os.getpid())

		tmp_pathway = tmp_files.add( 'pathway_pipeline.log' )
		PathwayPipeline(args.input_file, args.out_dir, tmp_pathway).submit(args.log_file)

	finally:
		print("Partie Finale ")


