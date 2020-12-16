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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse - Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
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

from frogsUtils import *
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class ParallelChimera(Cmd):
    """
    @summary: Removes PCR chimera by samples.
    """
    def __init__(self, in_fasta, in_abundance, out_fasta, out_abundance, out_summary, abundance_type, nb_cpus, log, debug, size_separator=None):
        """
        """
        size_separator_option = "" if size_separator is None else "--size-separator '" + size_separator + "' "
        debug_option = " --debug " if debug else ""
        Cmd.__init__( self,
                      'parallelChimera.py',
                      'Removes PCR chimera by samples.',
                      debug_option + size_separator_option + "--lenient-filter --nb-cpus " + str(nb_cpus) + " --sequences " + in_fasta + " --" + abundance_type + " " + in_abundance + " --non-chimera " + out_fasta + " --out-abundance " + out_abundance + " --summary " + out_summary + " --log-file " + log,
                      '--version' )


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_size_separator( in_fasta ):
    size_separator = ';size='
    FH_seq = FastaIO( in_fasta )
    record = FH_seq.next_seq()
    seq_idx = 1
    while seq_idx < 10 and size_separator is not None and record is not None:
        if not size_separator in record.id:
            size_separator = None
        elif size_separator is not None:
            try:
                real_id, size = record.id.rsplit(size_separator, 1)
                int(size)
            except:
                size_separator = None
        record = FH_seq.next_seq()
        seq_idx += 1
    FH_seq.close()
    return size_separator

def log_append_files( log_file, appended_files ):
    """
    @summary: Append content of several log files in one log file.
    @param log_file: [str] The log file where contents of others are appended.
    @param appended_files: [list] List of log files to append.
    """
    FH_log = Logger( log_file )
    FH_log.write( "\n" )
    for current_file in appended_files:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write( line )
        FH_input.close()
        FH_log.write( "\n" )
    FH_log.write( "\n" )
    FH_log.close()

def write_summary( summary_file, results_chimera ):
    """
    @summary: Writes the summary of results.
    @param summary_file: [str] The output file.
    @param results_chimera: [str] Path to the input chimera step summary.
    """
    # Get data
    detection_categories = ["Kept nb", "Kept abundance", "Removed nb", "Removed abundance", "Abundance of the most abundant removed", "Detected nb", "Detected abundance", "Abundance of the most abundant detected"]
    detection_data = list()
    remove_data = dict()

    # Parse results chimera
    in_remove_metrics = True
    in_detection_metrics = False
    section_first_line = True
    log_fh = open(results_chimera)
    for line in log_fh:
        line = line.strip()
        if line.startswith('##Metrics by sample'):
            remove_metrics = False
            in_detection_metrics = True
            section_first_line = True
        elif line.startswith('##Metrics global'):
            remove_metrics = True
            in_detection_metrics = False
            section_first_line = True
        elif line == "":
            in_detection_metrics = False
            in_remove_metrics = False
        else:
            if in_detection_metrics:
                if section_first_line:
                    line_fields = line[1:].split("\t")[1:]
                    detection_categories = line_fields
                    section_first_line = False
                else:
                    line_fields = line.split("\t")
                    detection_data.append({
                             'name': line_fields[0],
                             'data': list(map(int, line_fields[1:]))
                    })
            elif in_remove_metrics:
                if section_first_line:
                    line_fields = line[1:].split("\t")
                    remove_categories = [category.lower().replace(" ", "_") for category in line_fields]
                    section_first_line = False
                else:
                    for idx, val in enumerate(line.split("\t")):
                        remove_data[remove_categories[idx]] = int(val)
    log_fh.close()

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "remove_chimera_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###DETECTION_CATEGORIES###" in line:
            line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(detection_categories) )
        elif "###DETECTION_DATA###" in line:
            line = line.replace( "###DETECTION_DATA###", json.dumps(detection_data) )
        elif "###REMOVE_DATA###" in line:
            line = line.replace( "###REMOVE_DATA###", json.dumps(remove_data) )
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
        description='Removes PCR chimera.'
    )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-f', '--input-fasta', required=True, help='The cluster sequences (format: FASTA).' )
    group_exclusion_abundance = group_input.add_mutually_exclusive_group()
    group_exclusion_abundance.add_argument( '-b', '--input-biom', help='The abundance file for clusters by sample (format: BIOM).' )
    group_exclusion_abundance.add_argument( '-c', '--input-count', help='The counts file for clusters by sample (format: TSV).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-n', '--non-chimera', default='remove_chimera.fasta', help='sequences file without chimera (format: FASTA). [Default: %(default)s]')
    group_output.add_argument( '-a', '--out-abundance', default=None, help='Abundance file without chimera (format: BIOM or TSV). [Default: remove_chimera_abundance.biom or remove_chimera_abundance.tsv]')
    group_output.add_argument( '--summary', default="remove_chimera.html", help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.non_chimera)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

        tmp_chimera_summary = tmpFiles.add(os.path.basename(args.non_chimera) + "_summary.tsv")
        tmp_log  = tmpFiles.add(os.path.basename(args.non_chimera) + "_tmp.log")
        size_separator = get_size_separator( args.input_fasta )
        if args.input_count is None:
            if args.out_abundance == None:
                args.out_abundance = "remove_chimera_abundance.biom"
            ParallelChimera( args.input_fasta, args.input_biom, args.non_chimera, args.out_abundance, tmp_chimera_summary, "biom", args.nb_cpus, tmp_log, args.debug, size_separator ).submit( args.log_file )
        else:
            if args.out_abundance == None:
                args.out_abundance = "remove_chimera_abundance.tsv"
            ParallelChimera( args.input_fasta, args.input_count, args.non_chimera, args.out_abundance, tmp_chimera_summary, "count", args.nb_cpus, tmp_log, args.debug, size_separator ).submit( args.log_file )
        write_summary( args.summary, tmp_chimera_summary )
        
        # Append independant log files
        log_append_files( args.log_file, [tmp_log] )
        
    # Remove temporary files
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
