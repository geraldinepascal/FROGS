#!/usr/bin/env python2.7
#
# Copyright (C) 2015 INRA
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

__author__ = 'Maria Bernard - SIGENAE AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import argparse
from operator import itemgetter

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
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class SortFasta(Cmd):
    """
    @summary: Sort dereplicated sequences by decreasing abundance.
    """
    def __init__(self, in_fasta, out_fasta, size_separator=';size='):
        """
        @param in_fasta: [str] Path to unsorted file.
        @param out_fasta: [str] Path to file after sort.
        @param size_separator: [str] Each sequence in in_fasta is see as a pre-cluster. The number of sequences represented by the pre-cluster is stored in sequence ID.
               Sequence ID format : '<REAL_ID><size_separator><NB_SEQ>'. If this size separator is missing in ID, the number of sequences represented is 1.
        """
        Cmd.__init__( self,
                      'sortAbundancies.py',
                      'Sort pre-clusters by abundancies.',
                      "--size-separator '" + size_separator + "' --input-file " + in_fasta + ' --output-file ' + out_fasta,
                      '--version' )


class Swarm(Cmd):
    """
    @summary: Sequences clustering.
    @see: https://github.com/torognes/swarm
    """
    def __init__(self, in_fasta, out_swarms, out_log, distance, nb_cpus):
        """
        @param in_fasta: [str] Path to fasta file to process.
        @param out_swarms: [str] Path to swarm output file. It describes which reads compose each swarm.
        @param out_log: [str] Path to swarm log file.
        @param distance: [int] The 'param.distance'
        @param nb_cpus : [int] 'param.nb_cpus'.
        """
        Cmd.__init__( self,
                      'swarm',
                      'Clustering sequences.',
                      "--differences " + str(distance) + " --threads " + str(nb_cpus) + " --log " + out_log + " --output-file " + out_swarms + " " + in_fasta,
                      '--version' )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        try:
            p = Popen(self.program + ' ' + self.version_parameters, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = p.communicate()
            return stderr.split()[1]
        except:
            raise Exception( "Version cannot be retrieve for the software '" + self.program + "'." )


class Swarm2Biom(Cmd):
    """
    @summary: Converts swarm results in BIOM file.
    """
    def __init__(self, in_swarms, in_count, out_biom):
        """
        @param in_swarms: [str] Path to swarm output file. It describes which reads compose each swarm.
        @param in_count: [str] Path to the count file. It contains the count by sample for each representative sequence.
        @param out_biom: [str] Path to the output BIOM.
        """
        Cmd.__init__( self,
                      'swarm2biom.py',
                      'Converts swarm output to abundance file (format BIOM).',
                      "--clusters-file " + in_swarms + " --count-file " + in_count + " --output-file " + out_biom,
                      '--version' )


class ExtractSwarmsFasta(Cmd):
    """
    @summary: Extracts seeds sequences to produce the seeds fasta.
    """
    def __init__(self, in_fasta, in_swarms, out_seeds_file):
        """
        @param in_fasta: [str] Path to the input fasta file.
        @param in_swarms: [str] Path to swarm output file. It describes which reads compose each swarm.
        @param out_seeds_file: [str] Path to the output fasta file.
        """
        Cmd.__init__( self,
                      'extractSwarmsFasta.py',
                      'Extracts seeds sequences to produce the seeds fasta.',
                      '--input-fasta ' + in_fasta + ' --input-swarms ' + in_swarms + ' --output-fasta ' + out_seeds_file,
                      '--version' )


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
###
def resizeSeed(seed_in, seed_in_compo, seed_out):
    """
    @summary: add read abundance to seed sequence name
    @param seed_in : [str] Path to seed input fasta file
    @param seed_in_compo : [str] Path to seed input composition swarm file
    @param seed_out : [str] Path to seed output fasta file with abundance in name and sorted
    """
    dict_cluster_abond=dict()
    with open(seed_in_compo,"r") as f:
        for idx,line in enumerate(f.readlines()):
            if not line.startswith("#"):
                dict_cluster_abond["Cluster_"+str(idx+1)]=sum([ int(n.split("_")[-1]) for n in line.strip().split()])
    f.close()

    FH_input = FastaIO( seed_in )
    FH_out=FastaIO(seed_out , "w" )
    for record in FH_input:
        record.id += "_" + str(dict_cluster_abond[record.id])
        FH_out.write( record )
    FH_input.close()
    FH_out.close()

###
def agregate_composition(step1_compo , step2_compo, out_compo):
    """
    @summary: convert cluster composition in cluster in cluster composition in read (in case of two steps clustering)
    @param step1_compo : [str] Path to cluster1 composition in read (clustering step1)
    @param step2_compo : [str] Path to cluster2 composition in cluster1 (clustering step2) 
    @param out_composition : [str] Path to cluster2 composition in read
    """
    dict_cluster1_compo=dict()
    with open(step1_compo,"r") as f:
        for idx,line in enumerate(f.readlines()):
            dict_cluster1_compo["Cluster_"+str(idx+1)]=line.strip()
    f.close()

    FH_out=open(out_compo,"w")
    with open(step2_compo,"r") as f:
        for line in f.readlines():
            compo=" ".join([dict_cluster1_compo["_".join(n.split('_')[0:-1])] for n in line.strip().split(" ")])
            FH_out.write(compo+"\n")


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Single-linkage clustering on sequences.' )
    parser.add_argument( '-d', '--distance', type=int, default=2, help="Maximum distance between sequences in each aggregation step." )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used." )
    parser.add_argument( '-n', '--denoising', default=False, action='store_true',  help="denoise data by clustering read with distance=1 before perform real clustering" )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-f', '--input-fasta', required=True, help='The sequences file (format: fasta).' )
    group_input.add_argument( '-c', '--input-count', required=True, help="The count file for 'fasta-file' (format: count). It contains the count by sample for each sequence." )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-b', '--output-biom', default='swarms_abundance.biom', help='This output file will contain the abondance by sample for each cluster (format: BIOM).')
    group_output.add_argument( '--output-fasta', default='seeds.fasta', help='This output file will contain the seed sequence for each cluster (format: fasta).')
    group_output.add_argument( '--output-compo', default='swarms_composition.tsv', help='This output file will contain the composition of each cluster (format: TSV). One Line is a cluster ; each column is a sequence ID.')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )
    filename_woext = os.path.split(args.input_fasta)[1].split('.')[0]
    swarm_log = tmpFiles.add( filename_woext + '_swarm_log.txt' )
    sorted_fasta = tmpFiles.add( filename_woext + '_sorted.fasta' )
    final_sorted_fasta = sorted_fasta
    swarms_file = args.output_compo
    denoising_compo = None

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

        SortFasta( args.input_fasta, sorted_fasta ).submit( args.log_file )

        if args.denoising and args.distance > 1:
            # Denoising
            denoising_log = tmpFiles.add( filename_woext + '_denoising_log.txt' )
            denoising_compo = tmpFiles.add( filename_woext + '_denoising_composition.txt' )
            denoising_seeds = tmpFiles.add( filename_woext + '_denoising_seeds.fasta' )
            denoising_resized_seeds = tmpFiles.add( filename_woext + '_denoising_resizedSeeds.fasta' )
            swarms_file = tmpFiles.add( filename_woext + '_swarmD' + str(args.distance) + '_composition.txt' )
            final_sorted_fasta = tmpFiles.add( filename_woext + '_denoising_sortedSeeds.fasta' )
            Swarm( sorted_fasta, denoising_compo, denoising_log, 1 , args.nb_cpus ).submit( args.log_file )
            ExtractSwarmsFasta( sorted_fasta, denoising_compo, denoising_seeds ).submit( args.log_file )
            resizeSeed( denoising_seeds, denoising_compo, denoising_resized_seeds ) # add size to seeds name
            SortFasta( denoising_resized_seeds, final_sorted_fasta, "_" ).submit( args.log_file )

        Swarm( final_sorted_fasta, swarms_file, swarm_log, args.distance, args.nb_cpus ).submit( args.log_file )

        if args.denoising and args.distance > 1:
            # convert cluster composition in read composition ==> final swarm composition
            agregate_composition(denoising_compo, swarms_file, args.output_compo)

        Swarm2Biom( args.output_compo, args.input_count, args.output_biom ).submit( args.log_file )
        ExtractSwarmsFasta( final_sorted_fasta, swarms_file, args.output_fasta ).submit( args.log_file )

    # Remove temporary files
    finally:
        if not args.debug:
            tmpFiles.deleteAll()