#!/usr/bin/env python2.7
#
# Copyright (C) 2014 INRA
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

__author__ = 'Maria Bernard - SIGENAE'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import sys
import time
import json
import copy
import random
import argparse
from biom import *
from sequenceIO import *

def biom_fasta_update(biom_in, fasta_in, fasta_out, log_file):
    FH_in = FastaIO( fasta_in )
    FH_out = FastaIO( fasta_out, "w" )
    Biom = BiomIO.from_json( biom_in )
    seq_in=0
    seq_out=0

    for record in FH_in:
        seq_in += 1
        try:
            Biom.find_idx("observation",record.id)
        except ValueError:
            pass
        else:
            FH_out.write(record)
            seq_out += 1
    FH_in.close()
    FH_out.close()
    FH_log=open(log_file,"w")
    FH_log.write("Number of sequence in :" + str(seq_in)+"\n" )
    FH_log.write("Number of sequence out :" + str(seq_out) +"\n") 

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Update fasta file based on sequence in biom file' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-b', '--input-biom', required=True, help='Path to the BIOM input file.' )
    group_input.add_argument( '-f', '--input-fasta', required=True, help='Path to the Fasta input file.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='Path to the Fasta output file')
    group_output.add_argument( '-l', '--log-file', required=True, help='Path to the log output file')
    args = parser.parse_args()

    # Process
    biom_fasta_update(args.input_biom, args.input_fasta, args.output_file, args.log_file)