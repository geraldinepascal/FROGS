#!/usr/bin/env python3

__author__ = 'Frederic Escudie - Genotoul/MIAT & Maria Bernard - SIGENAE/GABI & Olivier Rué - Migale/MaIAGE'
__copyright__ = 'Copyright (C) 2025 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.1.0'
__email__ = 'frogs-support@inrae.fr'
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
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']
# THEME
THEME_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "static"))
if not os.path.exists(THEME_DIR):
    THEME_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(CURRENT_DIR)), "static"))

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
    def __init__(self, in_biom, out_newick, out_log, dist_method, linkage_method):
        """
        @param in_biom: [str] The processed BIOM path.
        @param out_newick: [str] The path to the output.
        @param out_log: [str] The path to the execution log.
        @param dist_method: [str] The distance method used.
        @param linkage_method: [str] The linkage method used.
        """
        self.exec_log = out_log
        Cmd.__init__( self,
                      'biomTools.py',
                      'Hierarchical classification on observation proportions.',
                      'hclassification --distance-method ' + dist_method + ' --linkage-method ' + linkage_method + ' --input-file ' + in_biom + ' --output-file ' + out_newick + ' > ' + out_log,
                      '--version' )

    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        excluded_exists = False

        # Parse execution log
        warning_lines = list()
        FH_exec_log = open( self.exec_log )
        for line in FH_exec_log:
            if line.strip() != "" and not line.startswith("#"):
                warning_lines.append(line.strip())
                if "xcluded samples" in line:
                    excluded_exists = True
        FH_exec_log.close()

        # Write warning (if at least one sample has been excluded)
        if excluded_exists:
            FH_log = Logger( log_file )
            FH_log.write( 'Warning:\n' )
            for line in warning_lines:
                FH_log.write( '\t' + line + '\n' )
            FH_log.close()
            
    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()            


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

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()                      


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_summary( summary_file, input_biom, depth_file, classif_file=None ):
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
    if classif_file is not None:
        FH_classif = open( classif_file )
        newick = FH_classif.readlines()[0].replace("\n", "")
        FH_classif.close()

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "cluster_asv_report_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    # Load shared JS
    with open(os.path.join(THEME_DIR, "js", "theme.js")) as f:
        theme_js = f.read()
    with open(os.path.join(THEME_DIR, "js", "utils.js")) as f:
        utils_js = f.read()
    # Load shared CSS
    with open(os.path.join(THEME_DIR, "css", "common.css")) as f:
        common_css = f.read()
    for line in FH_summary_tpl:
        if "###IMPORT_CSS###" in line:
            line = line.replace("###IMPORT_CSS###", f"<style type='text/css'>{common_css}</style>")
        elif "###IMPORT_JS_UTILS###" in line:
            # injection du JS inline
            line = line.replace("###IMPORT_JS_UTILS###", f"<script>\n{utils_js}</script>")
        elif "###IMPORT_JS_THEME###" in line:
            line = line.replace("###IMPORT_JS_THEME###", f"<script>\n{theme_js}</script>")
        elif "###CLUSTERS_SIZES###" in line:
            line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
        elif "###DATA_COUNTS###" in line:
            line = line.replace( "###DATA_COUNTS###", json.dumps(counts) )
        elif "###DATA_SAMPLE###" in line:
            line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
        elif "###NEWICK###" in line:
            if classif_file is not None:
                line = line.replace( "###NEWICK###", json.dumps(newick) )
            else:
                line = line.replace( "###NEWICK###", "null" )
        elif "###FROGS_VERSION###" in line:
            line = line.replace( "###FROGS_VERSION###", "\""+str(__version__)+"\"" )
        elif "###FROGS_TOOL###" in line:
            line = line.replace( "###FROGS_TOOL###", "\""+ os.path.basename(__file__)+"\"" )
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
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '--hierarchical-clustering', action='store_true', default=False, help="Perform Hierarchical classification on observation proportions. [Default: %(default)s]" )
    parser.add_argument( '--distance-method', type=str, default="braycurtis", help='Used distance method for classify (see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist). [Default: %(default)s]',
                         choices=["euclidean", "cityblock", "seuclidean", "sqeuclidean", "cosine", "correlation", "hamming", "jaccard", "chebyshev", "canberra", "braycurtis", "mahalanobis", "yule", "matching", "dice", "kulsinski", "rogerstanimoto", "russellrao", "sokalmichener", "sokalsneath", "wminkowski"] )
    parser.add_argument( '--linkage-method', type=str, default="average", help='Used linkage method for classify (see http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.linkage.html). [Default: %(default)s]',
                         choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"] )
    
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help='The BIOM file to process.' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--html', default='cluster_asv_report.html', help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument( '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)

    tmp_files = TmpFiles( os.path.split(args.html)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        if args.hierarchical_clustering:
            classif_file = tmp_files.add( "HClassif.newick" )
            classif_log = tmp_files.add( "HClassif_log.txt" )
            HClassification(args.input_biom, classif_file, classif_log, args.distance_method, args.linkage_method).submit( args.log_file )

        depth_file = tmp_files.add( "depths.tsv" )
        Depths(args.input_biom, depth_file).submit( args.log_file )
        if args.hierarchical_clustering:
            write_summary( args.html, args.input_biom, depth_file, classif_file )
        else:
            write_summary( args.html, args.input_biom, depth_file, None )

        
    # Remove temporary files
    finally:
        if not args.debug:
            tmp_files.deleteAll()
