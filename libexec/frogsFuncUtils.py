#!/usr/bin/env python3
#
# Copyright (C) 2022 INRA
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

__author__ = 'Vincent Darbot - GENPHYSE'
__copyright__ = 'Copyright (C) 2022 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse
import ete3 as ete

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import Biom, BiomIO
from frogsUtils import *
from frogsSequenceIO import *

####################################################################################################################
#
# Tree count
#
####################################################################################################################

def task_convert_fasta( args ):
	FH_input = FastaIO(args.input_fasta)
	FH_output = FastaIO(args.output_fasta, "wt")
	for record in FH_input:
		record.id = record.id
		record.description = None
		FH_output.write(record)
	FH_output.close()

def task_excluded_sequences( args ):
	"""
	@summary: Returns the excluded sequence, not insert into reference tree.
	@param fasta_file: [str] Path to the fasta file to process.
	@param tree_file: [str] Path to the tree file to process.
	@output: The file of no aligned sequence names.
	"""
	tree=ete.Tree(args.input_tree)
	all_leaves = list()
	for leaf in tree:
		all_leaves.append(leaf.name)

	FH_input = FastaIO(args.input_fasta)
	excluded = open(args.output_excluded, "wt")
	list_excluded = list()
	no_excluded = True 
	for record in FH_input:
		if record.id not in all_leaves:
			excluded.write(record.id+"\n")
			list_excluded.append(record.id)
			no_excluded = False
	FH_input.close()
	if no_excluded:
		excluded.write('#No excluded OTUs.\n')
	excluded.close()

####################################################################################################################
#
# Main
#
####################################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description='frogsFUncUtils is a program package designed for working with frogsfunc tools.' )
    parser.add_argument( "--version", action='version', version=__version__ )
    subparsers = parser.add_subparsers()

    # Convert fasta files
    parser_sampling = subparsers.add_parser('convert-fasta', help='Change fasta headers to be compatible with PICRUSt2.', usage='frogsFuncUtils.py convert-fasta [-h] -i INPUT_FILE -o OUTPUT_FILE')
    parser_sampling.add_argument( '-i', '--input-fasta', required=True, type=str, help='Fasta file processed.' )
    parser_sampling.add_argument( '-o', '--output-fasta', required=True, type=str, help='Output fasta file.' )   
    parser_sampling.set_defaults(func=task_convert_fasta)

    # Excluded sequences 
    parser_observationDepth = subparsers.add_parser('excluded-sequences', help='Returns the excluded sequence, not insert into reference tree.', usage='frogsFuncUtils.py convert-fasta [-h] -i INPUT_FILE -o OUTPUT_FILE')
    parser_observationDepth.add_argument( '-i', '--input-fasta', required=True, type=str, help='Input fasta file processed.' )
    parser_observationDepth.add_argument( '-t', '--input-tree', required=True, type=str, help='PICRUSt2 output tree with inserts sequences.' )
    parser_observationDepth.add_argument( '-e', '--output-excluded', required=True, type=str, help='Output file with sequences not inserts into reference tree.' )
    parser_observationDepth.set_defaults(func=task_excluded_sequences)
    # Parse parameters and call process
    args = parser.parse_args()
    args.func(args)
