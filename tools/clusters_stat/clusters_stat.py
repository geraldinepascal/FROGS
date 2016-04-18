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
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsBiom import BiomIO


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class HClassification(Cmd):
    """
    @summary: Hierarchical classification on observation proportions.
    """
    def __init__(self, in_biom, out_newick, dist_method, linkage_method):
        """
        @param in_biom: [str] The processed BIOM path.
        @param out_newick: [str] The path to the output.
        @param dist_method: [str] The distance method used.
        @param linkage_method: [str] The linkage method used.
        """
        Cmd.__init__( self,
                      'biomTools.py',
                      'Hierarchical classification on observation proportions.',
                      'hclassification --distance-method ' + dist_method + ' --linkage-method ' + linkage_method + ' --input-file ' + in_biom + ' --output-file ' + out_newick,
                      '--version' )


class Depths(Cmd):
    """
    @summary: Writes by abundance the number of clusters.
    """
    def __init__(self, in_biom, out_tsv):
        """
        @param in_biom: [str] The processed BIOM path.
        @param out_tsv: [str] The path of the output.
        """
        Cmd.__init__( self,
                      'biomTools.py',
                      'Writes by abundance the number of clusters.',
                      'obsdepth --input-file ' + in_biom + ' --output-file ' + out_tsv,
                      '--version' )


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_summary( summary_file, input_biom, depth_file, classif_file ):
    """
    @summary: Writes the summary of results.
    @param summary_file: [str] The output file.
    @param input_biom: [str] Path to the input BIOM.
    @param depth_file: [str] Path to biomTools obsdepth output.
    @param classif_file: [str] Path to biomTools hclassification output.
    """
    # Get size distribution data
    clusters_size = list()
    counts = list()
    FH_depth = open( depth_file )
    for line in FH_depth:
        if not line.startswith('#'):
            fields = line.strip().split()
            if fields[1] != "0":
                clusters_size.append( int(fields[0]) )
                counts.append( int(fields[1]) )
    FH_depth.close()

    # Get sample data
    biom = BiomIO.from_json( input_biom )
    samples_distrib = dict()
    for sample_name in biom.get_samples_names():
        shared_seq = 0
        shared_observations = 0
        own_seq = 0
        own_observations = 0
        for observation in biom.get_observations_by_sample(sample_name):
            obs_count_in_spl = biom.get_count( observation['id'], sample_name )
            if obs_count_in_spl != 0 and obs_count_in_spl == biom.get_observation_count(observation['id']):
                own_observations += 1
                own_seq += obs_count_in_spl
            else:
                shared_observations += 1
                shared_seq += obs_count_in_spl
        samples_distrib[sample_name] = {
            'shared_seq': shared_seq,
            'shared_observations': shared_observations,
            'own_seq': own_seq,
            'own_observations': own_observations
        }
    del biom

    # Get newick data
    FH_classif = open( classif_file )
    newick = FH_classif.readlines()[0].replace("\n", "")
    FH_classif.close()

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "clusters_stat_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###CLUSTERS_SIZES###" in line:
            line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
        elif "###DATA_COUNTS###" in line:
            line = line.replace( "###DATA_COUNTS###", json.dumps(counts) )
        elif "###DATA_SAMPLE###" in line:
            line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
        elif "###NEWICK###" in line:
            line = line.replace( "###NEWICK###", json.dumps(newick) )
        FH_summary_out.write( line )
    FH_summary_out.close()
    FH_summary_tpl.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Process several metrics on abundance from BIOM file.'
    )
    parser.add_argument( '--distance-method', type=str, default="braycurtis", help='Used distance method for classify (see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist).',
                         choices=["euclidean", "cityblock", "seuclidean", "sqeuclidean", "cosine", "correlation", "hamming", "jaccard", "chebyshev", "canberra", "braycurtis", "mahalanobis", "yule", "matching", "dice", "kulsinski", "rogerstanimoto", "russellrao", "sokalmichener", "sokalsneath", "wminkowski"] )
    parser.add_argument( '--linkage-method', type=str, default="average", help='Used linkage method for classify (see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.linkage.html).',
                         choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"] )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-i', '--input-biom', required=True, help='The BIOM file to process.' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-o', '--output-file', default='clusters_metrics.html', help='The output report.')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    tmp_files = TmpFiles( os.path.split(args.output_file)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

        classif_file = tmp_files.add( "HClassif.newick" )
        HClassification(args.input_biom, classif_file, args.distance_method, args.linkage_method).submit( args.log_file )

        depth_file = tmp_files.add( "depths.tsv" )
        Depths(args.input_biom, depth_file).submit( args.log_file )

        write_summary( args.output_file, args.input_biom, depth_file, classif_file )
    # Remove temporary files
    finally:
        if not args.debug:
            tmp_files.deleteAll()