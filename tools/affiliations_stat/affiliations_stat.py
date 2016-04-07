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
class Rarefaction(Cmd):
    """
    @summary: Writes by sample the rarefaction data.
    """
    def __init__(self, in_biom, tmp_files_manager, taxonomy_tag, rarefaction_levels):
        """
        @param in_biom: [str] The processed BIOM path.
        @param out_tsv: [str] The path of the output.
        @param taxonomy_tag: [str] The metadata title for the taxonomy in BIOM file.
        @param rarefaction_levels: [list] The taxonomy level(s) used to evaluate diversity.
        """
        # Step size management
        self.in_biom = in_biom
        step_size = self.get_step_size()
        # Out files management
        out_basename_pattern = "rarefaction_rank_##RANK##.tsv"
        out_files = list()
        for rank in rarefaction_levels:
            out_files.append( tmp_files_manager.add(out_basename_pattern.replace('##RANK##', str(rank))) )
        out_path_pattern = os.path.join( tmp_files_manager.tmp_dir, tmp_files_manager.prefix + "_" + out_basename_pattern )
        # Cmd
        Cmd.__init__( self,
                      'biomTools.py',
                      'Writes by sample the rarefaction data for rank(s) ' + ', '.join([str(lvl) for lvl in rarefaction_levels]) + '.',
                      'rarefaction --input-file ' + in_biom + ' --output-file-pattern ' + out_path_pattern + ' --taxonomy-key "' + taxonomy_tag + '" --step-size ' + str(step_size) + ' --ranks ' + ' '.join([str(lvl) for lvl in rarefaction_levels]),
                      '--version' )
        self.output_files = out_files

    def get_step_size(self, nb_step=35):
        """
        @summary: Returns the step size to obtain 'nb_step' steps or more in 3/4 of samples.
        @param nb_step: [int] The number of expected steps.
        @returns: [int] The step size.
        """
        counts = list()
        # Get the number of sequences by sample
        biom = BiomIO.from_json( self.in_biom )
        for sample_name in biom.get_samples_names():
            counts.append( biom.get_sample_count(sample_name) )
        del biom
        counts = sorted(counts)
        nb_samples = len(counts)
        # Finds the lower quartile number of sequences
        lower_quartile_idx = nb_samples/4
        nb_seq = counts[lower_quartile_idx]
        # If lower quartile sample is empty
        if nb_seq == 0:
            idx = 0
            while (lower_quartile_idx + idx) < nb_samples and counts[lower_quartile_idx + idx] == 0:
                nb_seq = counts[lower_quartile_idx + idx]
                idx += 1
        return int(nb_seq/nb_step)


class TaxonomyTree(Cmd):
    """
    @summary: Produces a taxonomy tree with counts by sample in extended newick format.
    """
    def __init__(self, in_biom, taxonomy_tag, out_tree, out_ids):
        """
        @param in_biom: [str] The processed BIOM path.
        @param taxonomy_tag: [str] The metadata title for the taxonomy in BIOM file.
        @param out_tree: [str] Path to the enewick output.
        @param out_ids: [str] Path to the IDs/samples output.
        """
        # Cmd
        Cmd.__init__( self,
                      'biomTools.py',
                      'Produces a taxonomy tree with counts by sample.',
                      'treeCount --input-file ' + in_biom + ' --taxonomy-key "' + taxonomy_tag + '" --output-enewick ' + out_tree + ' --output-samples ' + out_ids,
                      '--version' )


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_bootstrap_distrib( input_biom, bootstrap_tag, multiple_tag ):
    """
    @summary: Returns by taxonomic rank the count (seq and clstr) for the different bootstrap categories.
    @param input_biom: The path to the processed BIOM.
    @param bootstrap_tag: The metadata tag used in BIOM file to store the taxonomy bootstraps.
    @param multiple_tag: The metadata tag used in BIOM file to store the list of possible taxonomies.
    @returns: [dict] By taxonomic rank the count for the different bootstrap categories.
              Example:
                {
                    "Phylum": {
                        "80": { "clstr": 1, "seq":100 },
                        "90": {    "clstr": 2,    "seq":400 },
                        "100": { "clstr": 50, "seq":20000 },
                    },
                    "Genus":{
                        "80":{ "clstr": 1, "seq":100 },
                        "90":{ "clstr": 2, "seq":400 },
                        "100":{ "clstr": 50, "seq":20000 },
                    }
                }
    """
    bootstrap_results = dict()

    biom = BiomIO.from_json( input_biom )
    for observation in biom.get_observations():
        observation_metadata = observation['metadata']
        bootstrap = None
        if multiple_tag is not None:
            if observation_metadata.has_key(multiple_tag) and len(observation_metadata[multiple_tag]) > 0:
                bootstrap = observation_metadata[multiple_tag][0][bootstrap_tag]
        else:
            if observation_metadata.has_key(bootstrap_tag):
                bootstrap = observation_metadata[bootstrap_tag]
        if bootstrap is not None:
            for taxonomy_depth, rank_bootstrap in enumerate( bootstrap ):
                rank_bootstrap = rank_bootstrap * 100
                rank = args.taxonomic_ranks[taxonomy_depth]
                if not bootstrap_results.has_key(rank):
                    bootstrap_results[rank] = dict()
                if not bootstrap_results[rank].has_key(rank_bootstrap):
                    bootstrap_results[rank][rank_bootstrap] = {
                        "clstr": 0,
                        "seq": 0
                    }
                bootstrap_results[rank][rank_bootstrap]["clstr"] += 1
                bootstrap_results[rank][rank_bootstrap]["seq"] += biom.get_observation_count( observation['id'] )
    del biom
    return bootstrap_results

def get_alignment_distrib( input_biom, identity_tag, coverage_tag, multiple_tag ):
    """
    @summary: Returns by taxonomic rank the count (seq and clstr) for the different identity/coverage.
    @param input_biom: The path to the processed BIOM.
    @param identity_tag: The metadata tag used in BIOM file to store the alignment identity.
    @param coverage_tag: The metadata tag used in BIOM file to store the alignment query coverage.
    @param multiple_tag: The metadata tag used in BIOM file to store the list of possible taxonomies.
    @returns: [list] By taxonomic rank the count for the different identity/coverage.
              Example:
                [
                    [100, 100, { "clstr": 53, "seq": 20500 }],
                    [99, 100, { "clstr": 35, "seq": 18000 }],
                    [90, 95, { "clstr": 1, "seq": 10 }],
                ]
    """
    biom = BiomIO.from_json( input_biom )
    aln_results = list()
    aln_results_hash = dict()
    for observation in biom.get_observations():
        observation_metadata = observation['metadata']
        identity = None
        coverage = None
        if args.multiple_tag is not None:
            if observation_metadata.has_key(multiple_tag) and len(observation_metadata[multiple_tag]) > 0:
                identity = observation_metadata[multiple_tag][0][identity_tag]
                coverage = observation_metadata[multiple_tag][0][coverage_tag]
        else:
            if observation_metadata.has_key(identity_tag) and observation_metadata.has_key(coverage_tag):
                identity = observation_metadata[identity_tag]
                coverage = observation_metadata[coverage_tag]
        if identity is not None:
            if not aln_results_hash.has_key( identity ):
                aln_results_hash[identity] = dict()
            if not aln_results_hash[identity].has_key( coverage ):
                aln_results_hash[identity][coverage] = {
                    "clstr": 0,
                    "seq": 0
                }
            aln_results_hash[identity][coverage]["clstr"] += 1
            aln_results_hash[identity][coverage]["seq"] += biom.get_observation_count( observation['id'] )
    for ident in aln_results_hash.keys():
        for cover in aln_results_hash[ident].keys():
            aln_results.append([
                ident,
                cover,
                aln_results_hash[ident][cover]
            ])
    del biom
    return aln_results

def write_summary( summary_file, input_biom, tree_count_file, tree_ids_file, rarefaction_files, args ):
    """
    @summary: Writes the summary of results.
    @param summary_file: [str] The output file.
    @param input_biom: [str] Path to the input BIOM.
    @param tree_count_file: [str] Path to biomTools treeCount output.
    @param tree_ids_file: [str] Path to biomTools treeCount optional output.
    @param rarefaction_file: [str] Path to biomTools rarefaction output.
    @param args: The script arguments.
    """
    # Get taxonomy distribution
    FH_tree_count = open( tree_count_file )
    newick_tree = FH_tree_count.readline()
    FH_tree_count.close()
    ordered_samples_names = list()
    FH_tree_ids = open( tree_ids_file )
    for line in FH_tree_ids:
        id, sample_name = line.strip().split( "\t", 1 )
        ordered_samples_names.append( sample_name )
    FH_tree_ids.close()

    # Get bootstrap metrics
    bootstrap_results = None
    if args.bootstrap_tag is not None:
        bootstrap_results = get_bootstrap_distrib( input_biom, args.bootstrap_tag, args.multiple_tag )

    # Get alignment metrics
    aln_results = None
    if args.identity_tag is not None and args.coverage_tag is not None:
        aln_results = get_alignment_distrib( input_biom, args.identity_tag, args.coverage_tag, args.multiple_tag )

    # Get rarefaction data
    rarefaction_step_size = None
    rarefaction = None
    biom = BiomIO.from_json( input_biom )
    for rank_idx, current_file in enumerate(rarefaction_files):
        rank = args.rarefaction_ranks[rank_idx]
        FH_rarefaction = open( current_file )
        for line in FH_rarefaction:
            fields = map(str.strip, line.split("\t"))
            if line.startswith('#'):
                samples = fields[1:]
                if rarefaction is None:
                    rarefaction = dict()
                    for sample in samples:
                        rarefaction[sample] = dict()
                        rarefaction[sample]['nb_seq'] = biom.get_sample_count( sample )
                for sample in samples:
                    rarefaction[sample][rank] = list()
            else:
                if rarefaction_step_size is None:
                    rarefaction_step_size = int(fields[0])
                if not rarefaction[sample].has_key( rank ):
                    rarefaction[sample][rank] = list()
                for idx, sample in enumerate(samples):
                    if fields[idx+1] != "":
                        rarefaction[sample][rank].append( int(fields[idx+1]) )
        FH_rarefaction.close()
    del biom

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "affiliations_stat_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###TAXONOMIC_RANKS###" in line:
            line = line.replace( "###TAXONOMIC_RANKS###", json.dumps(args.taxonomic_ranks) )
        elif "###SAMPLES_NAMES###" in line:
            line = line.replace( "###SAMPLES_NAMES###", json.dumps(ordered_samples_names) )
        elif "###TREE_DISTRIBUTION###" in line:
            line = line.replace( "###TREE_DISTRIBUTION###", json.dumps(newick_tree) )
        elif "###DATA_RAREFACTION###" in line:
            line = line.replace( "###DATA_RAREFACTION###", json.dumps(rarefaction) )
        elif "###RAREFACTION_STEP_SIZE###" in line:
            line = line.replace( "###RAREFACTION_STEP_SIZE###", json.dumps(rarefaction_step_size) )
        elif "###RAREFACTION_RANKS###" in line:
            line = line.replace( "###RAREFACTION_RANKS###", json.dumps(args.rarefaction_ranks) )
        elif "###ALIGNMENT_SCORES###" in line:
            line = line.replace( "###ALIGNMENT_SCORES###", json.dumps(aln_results) )
        elif "###BOOTSTRAP_SCORES###" in line:
            line = line.replace( "###BOOTSTRAP_SCORES###", json.dumps(bootstrap_results) )
        FH_summary_out.write( line )
    FH_summary_out.close()
    FH_summary_tpl.close()

def process( args ):
    tmp_files = TmpFiles( os.path.split(args.output_file)[0] )

    try:
        # Add temp taxonomy if multiple and without consensus
        tmp_biom = args.input_biom
        used_taxonomy_tag = args.taxonomy_tag
        if args.multiple_tag is not None:
            used_taxonomy_tag = args.tax_consensus_tag
            if args.tax_consensus_tag is None:
                used_taxonomy_tag = "Used_taxonomy_FROGS-affi"
                tmp_biom = tmp_files.add( "tax.biom" )
                biom = BiomIO.from_json( args.input_biom )
                for observation in biom.get_observations():
                    metadata = observation["metadata"]
                    if len(metadata[args.multiple_tag]) > 0:
                        metadata[used_taxonomy_tag] = metadata[args.multiple_tag][0][args.taxonomy_tag]
                BiomIO.write( tmp_biom, biom )
                del biom

        # Rarefaction
        tax_depth = [args.taxonomic_ranks.index(rank) for rank in args.rarefaction_ranks]
        rarefaction_cmd = Rarefaction(tmp_biom, tmp_files, used_taxonomy_tag, tax_depth)
        rarefaction_cmd.submit( args.log_file )
        rarefaction_files = rarefaction_cmd.output_files

        # Taxonomy tree
        tree_count_file = tmp_files.add( "taxCount.enewick" )
        tree_ids_file = tmp_files.add( "taxCount_ids.tsv" )
        TaxonomyTree(tmp_biom, used_taxonomy_tag, tree_count_file, tree_ids_file).submit( args.log_file )

        # Writes summary
        write_summary( args.output_file, args.input_biom, tree_count_file, tree_ids_file, rarefaction_files, args )
    finally:
        if not args.debug:
            tmp_files.deleteAll()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description='Produces several metrics describing OTUs based on their taxonomies and the quality of the affiliations.')
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '--taxonomic-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels used in the metadata taxonomy.' )
    parser.add_argument( '--rarefaction-ranks', nargs='*', default=["Genus"], help='The ranks that will be evaluated in rarefaction.' )
    parser.add_argument( '--multiple-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the list of possible taxonomies. Use this parameter if the taxonomic affiliation has been processed by a software that adds several affiliation in the BIOM file (example: same score ambiguity).' )
    group_exclusion_taxonomy = parser.add_mutually_exclusive_group()
    group_exclusion_taxonomy.add_argument( '--tax-consensus-tag', type=str, help='The metadata tag used in BIOM file to store the consensus taxonomy. This parameter is used instead "--taxonomy-tag" when you have **************' )
    group_exclusion_taxonomy.add_argument( '--taxonomy-tag', type=str, help='The metadata tag used in BIOM file to store the taxonomy. Use this parameter if the taxonomic affiliation has been processed by a software that adds only one affiliation or if you does not have a metadata with the consensus taxonomy (see "--tax-consensus-tag").' )
    parser.add_argument( '--bootstrap-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the taxonomy bootstraps.' )
    parser.add_argument( '--identity-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the alignment identity.' )
    parser.add_argument( '--coverage-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the alignment observation coverage.' )
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-i', '--input-biom', required=True, help="The input biom file." )
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-o', '--output-file', default="affiliations_metrics.html", help="The output report." )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='The list of commands executed.' )
    args = parser.parse_args()
    prevent_shell_injections(args)

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")

    # Check parameters
    if args.multiple_tag is None and args.tax_consensus_tag is not None:
        raise Exception( "The parameter '--tax-consensus-tag' must be used only with the paameter '--multiple-tag'." )
    if args.taxonomy_tag is None and args.tax_consensus_tag is None:
        raise Exception( "The parameter '--taxonomy-tag' or the parameter '--tax-consensus-tag' must be set." )
    if (args.identity_tag is None and args.coverage_tag is not None) or (args.identity_tag is not None and args.coverage_tag is None):
        raise Exception( "The parameters '--identity-tag' and '--coverage-tag' must be setted together." )
    for current_rank in args.rarefaction_ranks:
        if current_rank not in args.taxonomic_ranks: raise Exception( "'" + current_rank + "' is not in valid taxonomic ranks : " + ", ".join(args.taxonomic_ranks) )
    biom = BiomIO.from_json( args.input_biom )
    if args.multiple_tag is None:
        for param in [args.taxonomy_tag, args.bootstrap_tag, args.identity_tag, args.coverage_tag]:
            if param is not None and not biom.has_observation_metadata( param ):
                raise Exception( "The metadata '" + param + "' does not exist in the BIOM file." )
    else:
        if args.tax_consensus_tag is not None and not biom.has_observation_metadata( args.tax_consensus_tag ):
            raise Exception( "The metadata '" + args.tax_consensus_tag + "' does not exist in the BIOM file." )
    del biom

    # Process
    process( args )