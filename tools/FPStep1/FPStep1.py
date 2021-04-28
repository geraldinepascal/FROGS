#!/usr/bin/env python3
#
# Copyright (C) 2018 INRA
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

__author__ = 'Maria Bernard INRA - SIGENAE AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2.1'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import json
import gzip
import argparse
import threading
import multiprocessing
import subprocess
from subprocess import Popen, PIPE

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
	def __init__(self, study_fasta, out_tree, ref_dir):
		"""
		@param study_fasta: [str] Path to input fasta file.
		@param out_tree: [str] Path to output resulting tree file.
		@param ref_dir: [str] Directory containing reference sequence files.
		"""
		if ref_dir != None:
			ref_dir = '--ref_dir '+ ref_dir
		else:
			ref_dir = ''

		Cmd.__init__(self,
		'place_seqs.py',
		'place OTU on reference tree.',
		'--study_fasta ' + study_fasta + ' --out_tree ' + out_tree + ref_dir + ' 2> stout.txt',
		'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').split()[1].strip()

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
	parser = argparse.ArgumentParser(description="place study sequences (i.e. OTUs) into a reference tree.")
	parser.add_argument('-v', '--version', action='version', version=__version__)
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-s', '--study_fasta', required=True,
	 help="Input fasta file of unaligned study sequences")
	group_input.add_argument('-r', '--ref_dir', default=None,
     help='Directory containing reference sequence files')
    # Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--out_tree', default='out.tree', help='Normalized sequences (format: FASTA).')
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')
	args = parser.parse_args()
	prevent_shell_injections(args)

	PlaceSeqs(args.study_fasta, args.out_tree, args.ref_dir).submit(args.log_file)




