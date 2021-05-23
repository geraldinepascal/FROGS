#!/usr/bin/env python3

__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os, sys
import argparse
import re

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
os.environ['PATH'] = CURRENT_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsBiom import BiomIO
from frogsSequenceIO import *

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def find_excluded(excluded_file):
	"""
	@summary: Returns the list of clusters not insert into reference tree from excluded.txt file.
	@param excluded_file: [str] Path to the excluded clusters file. 
	@note : Excluded clusters file must be one cluster ID per line:
	Cluster_1
	Cluster_4
	...
	"""
	excluded = []
	excluded_file = open(excluded_file,'r').readlines()
	if excluded_file is not None:
		for li in excluded_file:
			excluded.append(li.strip())
	return excluded

def remove_excluded_fasta( in_fasta, out_fasta, excluded_seqs):
	"""
	@summary: Returns the fasta file without sequences not insert into reference tree.
	"""
	FH_input = FastaIO( in_fasta )
	FH_output = FastaIO( out_fasta, "wt")
	for record in FH_input:
		real_id = record.id.split()[0]

		if real_id not in excluded_seqs:
			FH_output.write(record)
	FH_input.close()
	FH_output.close()

def remove_excluded_biom(input_biom, output_biom, excluded_seqs):
	"""
	@summary: Removes the specified list of observations.
	@param excluded_seqs: [list] The names of the observations to remove.
	"""
	biom = BiomIO.from_json( input_biom )
	biom.remove_observations( excluded_seqs )
	BiomIO.write( output_biom, biom )
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="Find OTUs closests references sequences in tree")
	parser.add_argument( '-v', '--version', action='version', version=__version__ )

	# Inputs
	group_input = parser.add_argument_group( 'Inputs' ) # Outputs
	group_input.add_argument( '-b', '--input_biom', required=True, help='Path to the abundance input file (format: biom).' )
	group_input.add_argument( '-i', '--input_fasta', required=True, help='Path to the sequences input file (format: fasta).' )
	group_input.add_argument( '-e', '--excluded_sequences', required=True, help='Path to the excluded observations IDs. Must be one ID per line.')
	# Outputs
	group_output = parser.add_argument_group( 'Outputs' ) # Outputs
	group_output.add_argument( '-o', '--output_biom', default='without_excluded_clusters.biom', help='Path to the Biom output file (format: biom).')
	group_output.add_argument( '-f', '--output_fasta', default='without_excluded_clusters.fasta', help='Path to the Fasta output file (format: fasta).')
	args = parser.parse_args()
	prevent_shell_injections(args)

	excluded = find_excluded(args.excluded_sequences)

	remove_excluded_fasta(args.input_fasta, args.output_fasta, excluded)
	remove_excluded_biom(args.input_biom, args.output_biom, excluded)