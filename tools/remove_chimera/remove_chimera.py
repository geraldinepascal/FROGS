#!/usr/bin/env python3

__author__ = 'Frederic Escudie - Genotoul/MIAT & Maria Bernard - SIGENAE/GABI'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

from decimal import DivisionByZero
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
from frogsBiom import Biom, BiomIO


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
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

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').strip()


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

def write_summary( summary_file, results_chimera, depth_file, biom_file):
    """
    @summary: Writes the summary of results.
    @param summary_file: [str] The output file.
    @param results_chimera: [str] Path to the input chimera step summary.
    """
    # Get data
    # detection_categories = ["Kept nb", "Kept abundance", "Removed nb", "Removed abundance", "Abundance of the most abundant removed", "Detected nb", "Detected abundance", "Abundance of the most abundant detected"]
    # detection_categories = ["Clusters kept", "Cluster abundance kept", "Chimeric clusters removed", "Chimeric abundance removed", "Abundance of the most abundant chimera removed", "Individual chimera detected", "Individual chimera abundance detected", "Abundance of the most abundant individual chimera detected"]
    detection_data = list()
    remove_data = dict()

    # Parse results chimera
    in_remove_metrics = True
    in_detection_metrics = False
    section_first_line = True
    log_fh = open(results_chimera)
    for line in log_fh:
        line = line.strip()
        if line.startswith('##Metrics global'):
            remove_metrics = True
            in_detection_metrics = False
            section_first_line = True
        elif line.startswith('##Metrics by sample'):
            remove_metrics = False
            in_detection_metrics = True
            section_first_line = True
        elif line == "":
            in_detection_metrics = False
            in_remove_metrics = False
        else:
            if in_detection_metrics:
                if section_first_line:
                    line_fields = line[1:].split("\t")[1:]
                    line_fields.insert(1,"%  Clusters kept")
                    line_fields.insert(3,"%  Cluster abundance kept")
                    detection_categories = line_fields
                    section_first_line = False
                else:
                    line_fields = line.split("\t")
                    # Calculate % of abundance kep on fly:
                    try:
                        line_fields.insert(2, float(100-(float(int(int(line_fields[3])*100)/\
                            int(int(line_fields[1])+int(line_fields[3]))))))
                    except ZeroDivisionError:
                        line_fields.insert(2, 'NA')
                    try:
                        line_fields.insert(4, float(100-(float(int(int(line_fields[5])*100)/\
                            int(int(line_fields[3])+int(line_fields[5]))))))
                    except ZeroDivisionError:
                        line_fields.insert(4, 'NA')     
                    detection_data.append({
                             'name': line_fields[0],
                             'data': line_fields[1:]
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
    biom = BiomIO.from_json( biom_file )
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
        elif "###DATA_SAMPLE###" in line:
            line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
        elif "###CLUSTERS_SIZES###" in line:
            line = line.replace( "###CLUSTERS_SIZES###", json.dumps(clusters_size) )
        elif "###DATA_COUNTS###" in line:
            line = line.replace( "###DATA_COUNTS###", json.dumps(counts) )
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
        description='Removes PCR chimera.'
    )
    parser.add_argument('--version', action='version', version=__version__ )
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )
    parser.add_argument('--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-fasta', required=True, help='The cluster sequences (format: FASTA).' )
    group_input.add_argument('--input-biom', required=True, help='The abundance file for clusters by sample (format: BIOM).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-fasta', default='remove_chimera.fasta', help='sequences file without chimera (format: FASTA). [Default: %(default)s]')
    group_output.add_argument('--output-biom', default='remove_chimera_abundance.biom', help='Abundance file without chimera (format: BIOM). [Default: %(default)s]')
    group_output.add_argument('--html', default="remove_chimera.html", help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.output_fasta)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

        tmp_chimera_summary = tmpFiles.add(os.path.basename(args.output_fasta) + "_summary.tsv")
        tmp_log  = tmpFiles.add(os.path.basename(args.output_fasta) + "_tmp.log")
        size_separator = get_size_separator( args.input_fasta )

        ParallelChimera( args.input_fasta, args.input_biom, args.output_fasta, args.output_biom, tmp_chimera_summary, "biom", args.nb_cpus, tmp_log, args.debug, size_separator ).submit( args.log_file )
        
        depth_file = tmpFiles.add( "depths.tsv" )
        Depths(args.output_biom, depth_file).submit( args.log_file )
        write_summary( args.html, tmp_chimera_summary, depth_file, args.output_biom )
        
        # Append independant log files
        log_append_files( args.log_file, [tmp_log] )
        
    # Remove temporary files
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
