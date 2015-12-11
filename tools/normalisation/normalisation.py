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

__author__ = 'Maria Bernard INRA - SIGENAE '
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.7.1'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'


import os
import sys
import json
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "bin"))
os.putenv('PATH', BIN_DIR + os.pathsep + os.getenv('PATH')) # $PATH
sys.path.insert(0, BIN_DIR) # $PYTHONPATH

from frogsUtils import *
from biom import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class BIOM_sampling(Cmd):
    """
    @summary: Random sampling in each sample.
    """
    def __init__(self, in_biom, out_biom, nb_read):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_biom: [str] Path to output BIOM file.
        @param nb_read : [int] Number of reads per sample.
        """
        Cmd.__init__( self,
                      "biomTools.py",
                      "Random sampling in each sample.",
                      "sampling --nb-sampled " + str(nb_read) + " --input-file " + in_biom + " --output-file " + out_biom,
                      "--version" )


class BIOM_FASTA_update(Cmd):
    """
    @summary: Converts BIOM file to TSV file.
    @note: taxonomyRDP seedID seedSequence blastSubject blastEvalue blastLength blastPercentCoverage blastPercentIdentity blastTaxonomy OTUname SommeCount sample_count
    """
    def __init__(self, in_biom, in_fasta, out_fasta, log):
        """
        @param in_biom: [str] Path to BIOM file.
        @param nb_read : [int] Number of reads per sample
        @param out_biom: [str] Path to output BIOM file.
        """
        Cmd.__init__( self,
                      'biomFastaUpdate.py',
                      'Update fasta file based on sequence in biom file',
                      "--input-biom " + in_biom + " --input-fasta " + in_fasta + " --output-file " + out_fasta + " --log " + log,
                      '--version' )


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_log(in_biom, out_biom, log):
    FH_log=open(log,"w")
    FH_log.write("#sample\tnb_otu_before\tnb_otu_after\n")
    initial_biom = BiomIO.from_json( in_biom )
    new_biom = BiomIO.from_json( out_biom )

    for sample_name in initial_biom.get_samples_names():
        nb_otu_before = len(initial_biom.get_sample_obs(sample_name))
        nb_otu_after = len(new_biom.get_sample_obs(sample_name))
        FH_log.write("Sample name: "+sample_name+"\n\tnb initials OTU: "+str(nb_otu_before)+"\n\tnb normalized OTU: "+str(nb_otu_after)+"\n")

    nb_initial_otu=len(initial_biom.rows)
    nb_new_otu=len(new_biom.rows)
    FH_log.write("Sample name: all samples\n\tnb initials OTU: "+str(nb_initial_otu)+"\n\tnb normalized OTU: "+str(nb_new_otu)+"\n")

    FH_log.close()

def summarise_results( summary_file, biom_subsample_log ):
    """
    @summary: Writes one summary of results from several logs.
    @param summary_file: [str] The output file.
    @param log_files: [list] The list of path to log files (one log file by sample).
    """
    # Get data
    categories = ["Nb OTU before normalisation" ,"Nb OTU after normalisation" ]
    series = list()
    get_sample_resuts(biom_subsample_log, series )
    histo = list()
    for i in xrange(0,len(series)) :
        if series[i]["name"]=="all samples":
            histo = series.pop(i)["data"]
            break
    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "normalisation_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###DATA_CATEGORIES###" in line:
            line = line.replace( "###DATA_CATEGORIES###", json.dumps(categories) )
        elif "###DATA_SERIES###" in line:
            line = line.replace( "###DATA_SERIES###", json.dumps(series) )
        elif "###HISTO_SERIES###" in line:
            line = line.replace( "###HISTO_SERIES###", json.dumps(histo) )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()

def get_sample_resuts( log_file, output_list ):
    """
    @summary: Returns the sample results (number of sequences after each filters).
    @param log_file: [str] Path to a log file.
    @return: [dict] The sample results. Format: {'name':SPL_NAME, 'data':[INI_NB_SEQ, NB_SEQ_AFTER_FILTER_1, NB_SEQ_AFTER_FILTER_2]}.
    """
    results = {
        'name': None,
        'data': list()
    }
    FH_input = open(log_file)
    for line in FH_input:
        if line.strip().startswith('Sample name: '):
            results['name'] = line.split(':')[1].strip()
        elif line.strip().startswith('nb initials OTU:'):
            results['data'].append( int(line.split(':')[1].strip()) )
        elif line.strip().startswith('nb normalized OTU: '):
            results['data'].append( int(line.split(':')[1].strip()) )
            output_list.append(results)
            results = {
                       'name': None,
                       'data': list()
                       }
    FH_input.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Normalisation in BIOM by random sampling.")
    parser.add_argument('-n', '--num-reads', type=int, required=True, help="Number of reads per sample after normalisation")
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-biom', required=True, help='Abundances file to normalize (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help='Sequences file to normalize (format: fasta).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-b', '--output-biom', default='abundance.biom', help='Normalized abundances (format: BIOM).')
    group_output.add_argument('-o', '--output-fasta', default='sequence.fasta', help='Normalized sequences (format: fasta).')
    group_output.add_argument('-s', '--summary-file', default='report.html', help='Summary of filters results (format: HTML).')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    tmp_files = TmpFiles( os.path.split(args.output_biom)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
        Logger.static_write(args.log_file,'Application start: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        
        Logger.static_write(args.log_file,'\n#Normalisation calculation\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        tmp_subsampling = tmp_files.add( 'tmp_biom_subsample.log' )
        BIOM_sampling(args.input_biom, args.output_biom, args.num_reads).submit(args.log_file)
        write_log(args.input_biom, args.output_biom, tmp_subsampling)
        Logger.static_write(args.log_file,'\tend: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
        tmp_fastaUpdate = tmp_files.add( 'tmp_fasta_update.log' )
        BIOM_FASTA_update(args.output_biom, args.input_fasta, args.output_fasta, tmp_fastaUpdate).submit(args.log_file)
        Logger.static_write(args.log_file,'\n#Summarise\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        summarise_results( args.summary_file, tmp_subsampling )
        Logger.static_write(args.log_file,'\tend: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
        Logger.static_write(args.log_file,'Application end: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
    # Remove temporary files
    finally:
        if not args.debug:
            tmp_files.deleteAll()
