#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import os
import sys
import argparse
import json
import re
from numpy import median

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

from frogsUtils import *
from frogsBiom import BiomIO
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES 
#
##################################################################################################################################################

class PlaceSeqs(Cmd):
	"""
	@summary: place study sequences (i.e. OTUs) into a reference tree
	"""
	def __init__(self, study_fasta, out_tree, placement_tool, ref_dir):
		"""
		@param study_fasta: [str] Path to input fasta file.
		@param out_tree: [str] Path to output resulting tree file.
		@param placement_tool: [str] Placement tool to use (epa-ng or sepp).
		@param ref_dir: [str] Directory containing reference sequence files.
		"""
		if ref_dir != None:
			ref_dir = ' --ref_dir '+ ref_dir
		else:
			ref_dir = ''

		Cmd.__init__(self,
		'place_seqs.py',
		'place OTU on reference tree.',
		'--study_fasta ' + study_fasta + ' --out_tree ' + out_tree + ' --placement_tool ' + placement_tool + ref_dir + ' 2> stout.txt',
		'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def convert_fasta(in_file, out_file):
	"""
	@summary: Change fasta headers to be compatible with picrust2
	"""
	FH_input = FastaIO(in_file)
	FH_output = FastaIO(out_file,"wt" )
	for record in FH_input:
		record.id = record.id
		record.description = None
		FH_output.write(record)
	FH_output.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="place study sequences (i.e. OTUs) into a reference tree.")
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-s', '--study_fasta', required=True,
	 help="Input fasta file of unaligned study sequences")
	group_input.add_argument('-r', '--ref_dir', default=None,
     help='Directory containing reference sequence files')
	group_input.add_argument('-t','--placement_tool', default='epa-ng',
	 help='Placement tool to use when placing sequences into reference tree. One of "epa-ng" or "sepp" must be input')
    # Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--out_tree', default='out.tree',
	 help='Normalized sequences (format: FASTA).')
	group_output.add_argument('-l', '--log-file', default=sys.stdout,
	 help='The list of commands executed.')
	args = parser.parse_args()
	prevent_shell_injections(args)

	tmpFiles=TmpFiles(os.path.split(args.out_tree)[0])

	try:

		fasta = tmpFiles.add('sout.fasta')
		convert_fasta(args.study_fasta,fasta)

		try:
			PlaceSeqs(args.study_fasta, args.out_tree, args.placement_tool, args.ref_dir).submit(args.log_file)

		except subprocess.CalledProcessError:
			print('\n\n#ERROR : epa-ng running out of memory. Please use placement tool sepp instead ( -t sepp )')
	finally:
		if not args.debug:
			tmp_files.deleteAll()
