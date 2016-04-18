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

__author__ = 'Katia Vidal - Team NED Toulouse AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.3.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'


import os
import sys
import json
import operator
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
from frogsSequenceIO import *


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class UpdateFasta(Cmd):
    """
    @summary: Updates fasta file based on sequence in biom file
    """
    def __init__(self, in_biom, in_fasta, out_fasta, log):
        """
        @param in_biom: [str] Path to BIOM file.
        @param nb_read : [int] Number of reads per sample
        @param out_biom: [str] Path to output BIOM file.
        """
        Cmd.__init__( self,
                      'biomFastaUpdate.py',
                      'Updates fasta file based on sequence in biom file.',
                      "--input-biom " + in_biom + " --input-fasta " + in_fasta + " --output-file " + out_fasta + " --log " + log,
                      '--version' )

class RemoveConta(Cmd):
    def __init__(self, contaminant_db, in_fasta, in_biom, cleaned_fasta, cleaned_biom, contaminated_fasta, log, nb_cpus ):
        """
        @param contaminant_db: [str] Path to the databank.
        @param in_fasta: [str] Path to the sequences file.
        @param in_biom: [str] Path to the abundance file.
        @param cleaned_fasta: [str] Path to the cleaned sequences file.
        @param cleaned_biom: [str] Path to the cleaned abundance file.
        @param contaminated_fasta: [str] Path to the contaminated sequences file.
        @param log: [str] Path to the log file.
        @param nb_cpus: [int] Number of usable CPUs.
        """
        Cmd.__init__( self,
                      "removeConta.py",
                      "Removes contaminant sequences.",
                      "--nb-cpus " + str(nb_cpus) + " --word-size 40 --min-identity 0.8 --min-coverage 0.8 --input-fasta " + in_fasta + " --contaminant-db " + contaminant_db + " --input-biom " + in_biom + " --clean-fasta " + cleaned_fasta + " --clean-biom " + cleaned_biom + " --conta-fasta " + contaminated_fasta + " --log-file " + log,
                      "--version")
        self.process_log = log
        self.contaminated_fasta = contaminated_fasta

    def write_conta_list(self, excluded_list):
        """
        @summary: Writes in output file the list of the sequences IDs tagged as contaminant.
        @excluded_list: [str] The path to the output file.
        """
        FH_contaminated_fasta = FastaIO( self.contaminated_fasta )
        FH_contaminated_list = open( excluded_list, "w" )
        for record in FH_contaminated_fasta:
            FH_contaminated_list.write( record.id + "\n" )
        FH_contaminated_fasta.close()
        FH_contaminated_list.close()

    def parser(self, log_file):
        """
        @summary: Parse the command results to add information in log_file.
        @log_file: [str] Path to the process log file.
        """
        FH_process_log = open(self.process_log, "r")
        for line in FH_process_log:
            if line.startswith("#Processed"):
                nb_processed = line.split(":")[1].strip()
            elif line.startswith("#Contaminated"):
                nb_contaminated = line.split(":")[1].strip()
        FH_process_log.close()
        FH_log = Logger(log_file)
        FH_log.write('Results:\n')
        FH_log.write('\tnumber of processed sequences: ' + str(nb_processed) + '\n')
        FH_log.write('\tnumber of removed contaminated sequences: ' + str(nb_contaminated) + '\n\n')
        FH_log.close()


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def excluded_obs_on_rdpBootstrap(input_biom, taxonomic_depth, min_bootstrap, excluded_file):
    """
    @summary: Writes the list of the observations with an insufficient bootstrap on the specified taxonomic rank.
    @param input_biom: [str] The path to the BIOM file to check.
    @param taxonomic_depth: [int] The taxonomic rank depth to check (example: 6 for Species in system "Domain, Phylum, Class, Order, Family, Genus, Species").
    @param min_bootstrap: [float] The observations with a value inferior to this threshold at the specified taxonomic depth are reported in the excluded file.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "w" )
    for observation in biom.get_observations():
        bootstrap = observation["metadata"]["rdp_bootstrap"]
        if issubclass(bootstrap.__class__, str):
            bootstrap = bootstrap.split(";")
        if bootstrap[taxonomic_depth] < min_bootstrap:
            FH_excluded_file.write( str(observation["id"]) + "\n" )
    FH_excluded_file.close()

def excluded_obs_on_samplePresence(input_biom, min_sample_presence, excluded_file):
    """
    @summary: Writes the list of the observations present in an insufficient number of samples.
    @param input_biom: [str] The path to the BIOM file to check.
    @param min_sample_presence: [int] The observations present in a number of samples inferior than this value are reported in the excluded file.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "w" )
    for observation_name in biom.get_observations_names():
        nb_samples = sum(1 for x in biom.get_samples_by_observation(observation_name))
        if nb_samples < min_sample_presence:
            FH_excluded_file.write( observation_name + "\n" )
    FH_excluded_file.close()

def excluded_obs_on_abundance(input_biom, min_abundance, excluded_file):
    """
    @summary: Writes the list of the observations with an insufficient abundance.
    @param input_biom: [str] The path to the BIOM file to check.
    @param min_abundance: [int/float] The observations with an abundance inferior than this value are reported in the excluded file.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "w" )
    min_nb_seq = min_abundance
    if type(min_abundance) == float:
        min_nb_seq = biom.get_total_count() * min_abundance
    for idx, count_by_sample in enumerate(biom.to_count()):
        observation = biom.rows[idx]
        abundance = sum(count_by_sample)
        if abundance < min_nb_seq:
            FH_excluded_file.write( str(observation["id"]) + "\n" )
    FH_excluded_file.close()

def excluded_obs_on_blastMetrics( input_biom, tag, cmp_operator, threshold, excluded_file ):
    """
    @summary: Writes the list of the observations with no affiliations with sufficient blast value.
    @param input_biom: [str] The path to the BIOM file to check.
    @param tag: [str] The metadata checked.
    @param cmp_operator: [str] The operator use in comparison (tag_value ">=" thresold or tag_value "<=" thresold ).
    @param threshold: [float] The limit for the tag value.
    @param excluded_file: [str] The path to the output file.
    """
    valid_operators = {
        ">=": operator.__ge__,
        "<=": operator.__le__
    }
    cmp_func = valid_operators[cmp_operator]
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "w" )
    for observation in biom.get_observations():
        alignments = observation["metadata"]["blast_affiliations"]
        is_discarded = True
        for current_alignment in alignments:
            if cmp_func(float(current_alignment[tag]), threshold):
                is_discarded = False
        if is_discarded:
            FH_excluded_file.write( str(observation["id"]) + "\n" )
    FH_excluded_file.close()

def excluded_obs_on_nBiggest( input_biom, nb_selected, excluded_file ):
    """
    @summary: Writes the list of all the observations without the n most abundant.
    @param input_biom: [str] The path to the BIOM file.
    @param threshold: [float] The number of the most abundant observations that will not be written in the excluded list.
    @param excluded_file: [str] The path to the output file.
    """
    biom = BiomIO.from_json( input_biom )
    FH_excluded_file = open( excluded_file, "w" )
    sorted_obs_counts = sorted( biom.get_observations_counts(), key=lambda observation: observation[1], reverse=True )
    for observation_name, observation_count in sorted_obs_counts[nb_selected:]:
        FH_excluded_file.write( observation_name + "\n" )
    FH_excluded_file.close()

def remove_observations( removed_observations, input_biom, output_biom ):
    """
    @summary: Removes the specified list of observations.
    @param removed_observations: [list] The names of the observations to remove.
    @param input_biom: [str] The path to the input BIOM.
    @param output_biom: [str] The path to the output BIOM.
    """
    biom = BiomIO.from_json( input_biom )
    biom.remove_observations( removed_observations )
    BiomIO.write( output_biom, biom )

def uniq_from_files_lists( in_files ):
    """
    @summary: Returns an list without duplicated elements from several list files.
    @param in_files: [list] The list of files paths. Each file contains a list.
    @return: [list] The list without duplicated elements.
    """
    uniq = dict()
    for current_file in in_files:
        FH_current_file = open( current_file )
        for line in FH_current_file:
            uniq[line.strip()] = 1
        FH_current_file.close()
    return uniq.keys()

def write_exclusion( discards, excluded_file ):
    """
    @summary: Writes the list of observations removed by each filter.
    @param discards: [dict] By filter the path of the file that contains the list of the removed observations.
    @param excluded_file: [str] The path to the output file.
    """
    FH_excluded = open( excluded_file, "w" )
    list_FH_discards = list()

    # Header
    header_line_fields = list()
    for filter in discards:
        header_line_fields.append( filter )
        list_FH_discards.append( open(discards[filter]) )
    FH_excluded.write( "#" + "\t".join(header_line_fields)  + "\n" )

    # Excluded
    nb_eof = 0
    while nb_eof < len(list_FH_discards):
        discards_line_fields = list()
        for FH_idx, FH_curent_filter in enumerate(list_FH_discards): # For each filter
            observation = ""
            if FH_curent_filter is not None: # Process next line if the discard file is not closed
                observation = FH_curent_filter.readline()
                if observation == "":
                    FH_curent_filter.close()
                    list_FH_discards[FH_idx] = None
                    nb_eof += 1
            discards_line_fields.append( observation.strip() )
        if nb_eof < len(list_FH_discards):
            FH_excluded.write( "\t".join(discards_line_fields)  + "\n" )
    FH_excluded.close()

def write_summary( summary_file, input_biom, output_biom, discards ):
    """
    @summary: Writes the process summary.
    @param summary_file: [str] The path to the output file.
    @param input_biom: [str] The path to the BIOM before program execution.
    @param output_biom: [str] The path to the BIOM after program execution.
    @param discards: [dict] By filter the path of the file that contains the list of the removed observations.
    """
    global_results = {
        'nb_clstr_kept': 0,
        'nb_clstr_ini': 0,
        'nb_seq_kept': 0,
        'nb_seq_ini': 0
    }
    samples_results = dict()
    filters_results = dict()

    # Global before filters
    in_biom = BiomIO.from_json( input_biom )
    for observation_name in in_biom.get_observations_names():
        global_results['nb_clstr_ini'] += 1
        global_results['nb_seq_ini'] += in_biom.get_observation_count( observation_name )
    for sample_name in in_biom.get_samples_names():
        samples_results[sample_name] = {
            'initial': sum( 1 for x in in_biom.get_observations_by_sample(sample_name) ),
            'filtered': dict(),
            'kept': 0
        }

    # By sample and by filters
    filters_intersections = dict()
    for filter in discards.keys():
        FH_filter = open( discards[filter] )
        for line in FH_filter:
            observation_name = line.strip()
            if not filters_intersections.has_key( observation_name ):
                filters_intersections[observation_name] = dict()
            filters_intersections[observation_name][filter] = 1
        FH_filter.close()
    for observation_name in filters_intersections.keys():
        # Removed intersection
        intersections_key = "--@@--".join(sorted( filters_intersections[observation_name].keys() ))
        if not filters_results.has_key( intersections_key ):
            filters_results[intersections_key] = {
                'filters': filters_intersections[observation_name].keys(),
                'count': 0
            }
        filters_results[intersections_key]['count'] += 1

        # Filters by samples
        for sample in in_biom.get_samples_by_observation(observation_name):
            for filter in filters_intersections[observation_name]:
                if not samples_results[sample['id']]['filtered'].has_key(filter):
                    samples_results[sample['id']]['filtered'][filter] = 0
                samples_results[sample['id']]['filtered'][filter] += 1
    del in_biom

    # Global after filters
    out_biom = BiomIO.from_json( output_biom )
    for observation_name in out_biom.get_observations_names():
        global_results['nb_clstr_kept'] += 1
        global_results['nb_seq_kept'] += out_biom.get_observation_count( observation_name )
    for sample_name in out_biom.get_samples_names():
        samples_results[sample_name]['kept'] = sum( 1 for x in out_biom.get_observations_by_sample(sample_name) )
    del out_biom

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "filters_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###PORCESSED_FILTERS###" in line:
            line = line.replace( "###PORCESSED_FILTERS###", json.dumps([filter for filter in discards]) )
        elif "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
        elif "###FILTERS_RESULTS###" in line:
            line = line.replace( "###FILTERS_RESULTS###", json.dumps(filters_results.values()) )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()

def ratioParameter( arg_value ):
    """
    @summary: Argparse type for ratio (float between 0 and 1).
    """
    float_arg_value = None
    try:
        float_arg_value = float(arg_value)
        if float_arg_value < 0.0 or float_arg_value > 1.0:
            raise argparse.ArgumentTypeError("must be between 0.0 and 1.0.")
    except:
        raise argparse.ArgumentTypeError("must be between 0.0 and 1.0.")
    return float_arg_value

def minAbundParameter( arg_value ):
    """
    @summary: Argparse type for min-abundance parameter.
    """
    str_value = arg_value.replace(",",".") 
    cleaned_value = float(str_value) if "." in str_value or "e-" in str_value else int(str_value)
    return cleaned_value

class BootstrapParameter(argparse.Action):
    """
    @summary : Argparse parameter for min-rdp-bootstrap parameter.
    """
    def __call__(self, parser, namespace, value, option_string=None):
        # Set parser
        output = getattr(namespace, self.dest)
        if output is None:
            output = {
                "rank": None,
                "value": None
            }
        if len(value.split(":")) != 2:
            raise argparse.ArgumentTypeError("The parameter '--min-rdp-bootstrap' must be in format 'TAXONOMIC_LEVEL:MIN_BOOTSTRAP'.")
        output["rank"] = value.split(":")[0]
        output["value"] = value.split(":")[1]
        try:
            output["value"] = ratioParameter(output["value"])
        except:
            raise argparse.ArgumentTypeError("The value for the MIN_BOOTSTRAP in parameter '--min-rdp-bootstrap' must be between 0.0 and 1.0.")
        setattr(namespace, self.dest, output)

def process( args ):
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )

    try:
        discards = dict() # by filter the discard file path

        if args.min_sample_presence is not None:
            label = "Present in minus of " + str(args.min_sample_presence) + " samples"
            discards[label] = tmpFiles.add( "min_sample_presence" )
            excluded_obs_on_samplePresence( args.input_biom, args.min_sample_presence, discards[label] )

        if args.min_abundance is not None:
            label = "Abundance < " + str(args.min_abundance)
            discards[label] = tmpFiles.add( "min_abundance" )
            excluded_obs_on_abundance( args.input_biom, args.min_abundance, discards[label] )

        if args.min_rdp_bootstrap is not None:
            label = "RDP bootstrap for " + args.min_rdp_bootstrap["rank"] + " < " + str(args.min_rdp_bootstrap["value"])
            discards[label] = tmpFiles.add( "min_rdp_bootstrap" )
            excluded_obs_on_rdpBootstrap( args.input_biom, args.rdp_taxonomy_ranks.index(args.min_rdp_bootstrap["rank"]), args.min_rdp_bootstrap["value"], discards[label] )

        if args.min_blast_length is not None:
            label = "All blast length < " + str(args.min_blast_length)
            discards[label] = tmpFiles.add( "min_blast_length" )
            excluded_obs_on_blastMetrics( args.input_biom, "aln_length", ">=", args.min_blast_length, discards[label] )

        if args.max_blast_evalue is not None:
            label = "All blast evalue > " + str(args.max_blast_evalue)
            discards[label] = tmpFiles.add( "max_blast_evalue" )
            excluded_obs_on_blastMetrics( args.input_biom, "evalue", "<=", args.max_blast_evalue, discards[label] )

        if args.min_blast_identity is not None:
            label = "All blast identity < " + str(args.min_blast_identity)
            discards[label] = tmpFiles.add( "min_blast_identity" )
            excluded_obs_on_blastMetrics( args.input_biom, "perc_identity", ">=", 100*args.min_blast_identity, discards[label] )

        if args.min_blast_coverage is not None:
            label = "All blast coverage < " + str(args.min_blast_coverage)
            discards[label] = tmpFiles.add( "min_blast_coverage" )
            excluded_obs_on_blastMetrics( args.input_biom, "perc_query_coverage", ">=", 100*args.min_blast_coverage, discards[label] )

        if args.contaminant is not None:
            label = "Present in databank of contaminants"
            discards[label] = tmpFiles.add( "contaminant" )
            cleaned_seq = tmpFiles.add( "cleaned_sequences.fasta" )
            cleaned_biom = tmpFiles.add( "cleaned_abundance.biom" )
            contaminated_seq = tmpFiles.add( "contaminated_sequences.fasta" )
            remove_log = tmpFiles.add( "clean.log" )
            cmd = RemoveConta(args.contaminant, args.input_fasta, args.input_biom, cleaned_seq, cleaned_biom, contaminated_seq, remove_log, args.nb_cpus)
            cmd.submit( args.log_file )
            cmd.write_conta_list( discards[label] )

        # Removes observations
        clusters_discarded = uniq_from_files_lists( [discards[filter] for filter in discards] )
        remove_observations( clusters_discarded, args.input_biom, args.output_biom )

        # Selects the N most abundant
        if args.nb_biggest_otu:
            label = "Not in the " + str(args.nb_biggest_otu) + " biggest"
            discards[label] = tmpFiles.add( "nb_biggest_otu" )
            excluded_obs_on_nBiggest( args.output_biom, args.nb_biggest_otu, discards[label] )
            # Removes observations
            remove_observations( uniq_from_files_lists([discards[label]]), args.output_biom, args.output_biom )

        # Writes outputs
        update_fasta_log = tmpFiles.add( "update_fasta_log.txt" )
        UpdateFasta( args.output_biom, args.input_fasta, args.output_fasta, update_fasta_log ).submit( args.log_file )
        write_exclusion( discards, args.excluded )
        write_summary( args.summary, args.input_biom, args.output_biom, discards )

    finally:
        tmpFiles.deleteAll()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description='Filters an abundance file')
    parser.add_argument('-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used.")
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    #     Filters
    group_filter = parser.add_argument_group( 'Filters' )
    group_filter.add_argument( '--nb-biggest-otu', type=int, default=None, required=False, help="Number of most abundant OTUs you want to keep.") 
    group_filter.add_argument( '-s', '--min-sample-presence', type=int, help="Keep OTU present in at least this number of samples.") 
    group_filter.add_argument( '-a', '--min-abundance', type=minAbundParameter, default=None, required=False, help="Minimum percentage/number of sequences, comparing to the total number of sequences, of an OTU (between 0 and 1 if percentage desired)." )  
    group_filter.add_argument( '-b', '--min-rdp-bootstrap', type=str, action=BootstrapParameter, metavar=("TAXONOMIC_LEVEL:MIN_BOOTSTRAP"), help="The minimal RDP bootstrap must be superior to this value (between 0 and 1)." )
    group_filter.add_argument( '-t', '--rdp-taxonomy-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels present in the reference databank.' )
    group_filter.add_argument( '-i', '--min-blast-identity', type=ratioParameter, help="The number corresponding to the blast percentage identity (between 0 and 1)." )
    group_filter.add_argument( '-c', '--min-blast-coverage', type=ratioParameter, help="The number corresponding to the blast percentage coverage (between 0 and 1)." )
    group_filter.add_argument( '-e', '--max-blast-evalue', type=float, help="The number corresponding to the blast e value (between 0 and 1).")
    group_filter.add_argument( '-l', '--min-blast-length', type=int, default=None, required=False, help="The number corresponding to the blast length." )
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help="The input biom file.")  #biomTest.txt
    group_input.add_argument('--input-fasta', required=True, help="The input fasta file.")  #test.fasta
    group_input.add_argument('--contaminant', default=None, help="Use this databank to filter sequence before affiliation.")
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default="filtered.biom", help="The Biom file output.")
    group_output.add_argument('--output-fasta', default="filtered.fasta", help="The fasta output file.")
    group_output.add_argument('--summary', default="summary.html", help="The HTML file containing the graphs.")
    group_output.add_argument('--excluded', default="excluded.tsv", help="The file that summarizes all the clusters discarded.")
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")

    if args.nb_biggest_otu is None and args.min_sample_presence is None and args.min_abundance is None and args.min_rdp_bootstrap is None and args.min_blast_identity is None and args.min_blast_coverage is None and args.max_blast_evalue is None and args.min_blast_length is None:
        raise argparse.ArgumentTypeError( "At least one filter must be set to run " + os.path.basename(sys.argv[0]) )
    if args.min_abundance <= 0 or (type(args.min_abundance) == float and args.min_abundance >= 1.0 ) :
        raise argparse.ArgumentTypeError( "If filtering on abundance, you must indicate a positiv threshold and if percentage abundance threshold must be smaller than 1.0. " )
    in_biom = BiomIO.from_json( args.input_biom )
    if args.min_rdp_bootstrap is not None:
        if not in_biom.has_observation_metadata("rdp_bootstrap"):
            raise argparse.ArgumentTypeError( "The BIOM input does not contain the metadata 'rdp_bootstrap'. You cannot use the parameter '--min-rdp-bootstrap' on this file." )
        if not args.min_rdp_bootstrap["rank"] in args.rdp_taxonomy_ranks:
            raise argparse.ArgumentTypeError( "The taxonomic rank choosen in '--min-rdp-bootstrap' must be in '--rdp-taxonomy-ranks' (" + ";".join(args.rdp_taxonomy_ranks) + ")." )
    if not in_biom.has_observation_metadata("blast_affiliations"):
        if args.min_blast_length is not None:
            raise argparse.ArgumentTypeError( "The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-length' on this file." )
        if args.max_blast_evalue is not None:
            raise argparse.ArgumentTypeError( "The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-evalue' on this file." )
        if args.min_blast_identity is not None:
            raise argparse.ArgumentTypeError( "The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-identity' on this file." )
        if args.min_blast_coverage is not None:
            raise argparse.ArgumentTypeError( "The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-coverage' on this file." )
    del in_biom

    # Process
    process( args )
