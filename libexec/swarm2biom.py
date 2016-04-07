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
import time
import json
import copy
import random
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import Biom, BiomIO


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def to_biom( clusters_file, count_file, output_biom, size_separator ):
    """
    @summary : Write a biom file from swarm results.
    @param clusters_file : [str] path to the '.clstr' file.
    @param count_file : [str] path to the count file. It contains the count of
                         sequences by sample of each preclusters.
                         Line format : "Precluster_id    nb_in_sampleA    nb_in_sampleB"
    @param output_biom : [str] path to the output file.
    @param size_separator : [str] the pre-cluster abundance separator.
    """
    biom = Biom( generated_by='swarm', matrix_type="sparse" )

    # Preclusters count by sample
    preclusters_count = dict()
    count_fh = open( count_file )
    samples = count_fh.readline().strip().split()[1:]
    for line in count_fh:
        line_fields = line.strip().split()
        count_by_sample = {}
        for idx, val in enumerate(line_fields[1:]):
            if val > 0:
                count_by_sample[samples[idx]] = int(val)
        preclusters_count[line_fields[0]] = count_by_sample
    count_fh.close()

    # Add samples
    for sample_name in samples:
        biom.add_sample( sample_name )

    # Process count
    cluster_idx = 1
    clusters_fh = open( clusters_file )
    for line in clusters_fh:
        cluster_name = "Cluster_" + str(cluster_idx)
        cluster_count = {key:0 for key in samples}
        line_fields = line.strip().split()
        # Retrieve count by sample
        for seq_id in line_fields:
            real_seq_id = seq_id.rsplit(size_separator, 1)[0]
            for preclust_sample in preclusters_count[real_seq_id]:
                cluster_count[preclust_sample] += preclusters_count[real_seq_id][preclust_sample]
            preclusters_count[real_seq_id] = None
        # Add cluster on biom
        biom.add_observation( cluster_name, {'seed_id':line_fields[0].rsplit(size_separator, 1)[0]} )
        for sample_name in samples:
            if cluster_count[sample_name] > 0:
                biom.add_count( cluster_name, sample_name, cluster_count[sample_name] )
        # Next cluster
        cluster_idx += 1

    # Write
    BiomIO.write( output_biom, biom )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Convert Swarm results to Biom file.' )
    parser.add_argument( '-s', '--size-separator', default='_', action='store_true', help="The pre-cluster abundance separator." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '--clusters-file', required=True, help='The cluster file from swarm (format : TSV).' )
    group_input.add_argument( '--count-file', required=True, help='Path to the count file (format : TSV). It contains the count of sequences by sample of each preclusters.' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', default='swarm_abundance.biom', help='Path to the output file (format : BIOM).')
    args = parser.parse_args()

    # Process
    to_biom( args.clusters_file, args.count_file, args.output_file, args.size_separator )
