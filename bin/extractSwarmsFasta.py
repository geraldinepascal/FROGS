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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import gzip
import argparse
from sequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def filter_seq( input_fasta, clusters_file, cluster_fasta ):
    """
     @summary: Write a renamed_fasta where the representative sequences ID will be replaced by the ID of the cluster.
      @param input_fasta : [str] path to the fasta to process.
      @param clusters_file : [str] path to the '.clstr'.
      @param renamed_fasta : [str] path to the fasta after process.
    """
    cluster_representative = dict()
    # Retrieve representatives sequences
    cluster_idx = 1
    clusters_fh = open( clusters_file )
    for line in clusters_fh:
        line = line.strip()
        representative_id = line.split()[0]
        cluster_representative[representative_id] = "Cluster_" + str(cluster_idx)
        cluster_idx += 1
    clusters_fh.close()

    # Filter sequences
    FH_in = FastaIO( input_fasta )
    FH_out = FastaIO( cluster_fasta, "w" )
    for record in FH_in:
        if cluster_representative.has_key(record.id):
            record.id = cluster_representative[record.id]
            FH_out.write( record )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Extracts seeds sequences to produce the seeds fasta.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-f', '--input-fasta', required=True, help='This file cotains sequences used in clustering  (format : Fasta).' )
    group_input.add_argument( '-s', '--input-swarms', required=True, help='This file contains the composition of each cluster (format : TSV). One Line is a cluster ; each column is a sequence ID.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-fasta', required=True, help='This file will contain the seed sequence for each cluster (format : Fasta).')
    args = parser.parse_args()

    # Process
    filter_seq( args.input_fasta, args.input_swarms, args.output_fasta )