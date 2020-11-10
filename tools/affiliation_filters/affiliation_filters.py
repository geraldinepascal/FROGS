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

__author__ = 'Maria Bernard - SIGENAE - Jouy-en-Josas'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'


import os
import sys
import copy
import json
import operator
import argparse
import re

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
APP_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "app"))
os.environ['PATH'] = APP_DIR + os.pathsep + os.environ['PATH']

# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

from frogsUtils import *
from frogsBiom import *
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

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

class Biom_to_tsv(Cmd):
    """
    @summary: Convert a biom abundance file into TSV files (abundance and multihit)
    """
    def __init__(self, in_biom, out_tsv, out_multihit):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv : [int] Path to abundance TSV file 
        @param out_multihit: [str] Path to multi affiliation TSV file
        """
        Cmd.__init__( self,
                      'biom_to_tsv.py',
                      'Convert a biom abundance file into TSV files (abundance and multihit)',
                      "--input-biom " + in_biom + " --output-tsv " + out_tsv + " --output-multi-affi " + out_multihit,
                      '--version' )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def ratioParameter( arg_value ):
    """
    @summary: Argparse type for ratio (float between 0 and 1).
    """
    float_arg_value = None
    try:
        float_arg_value = float(arg_value)
        if float_arg_value < 0.0 or float_arg_value > 1.0:
            raise argparse.ArgumentTypeError("\n\n#ERROR : must be between 0.0 and 1.0.\n")
    except:
        raise argparse.ArgumentTypeError("\n\n#ERROR : must be between 0.0 and 1.0.\n")
    return float_arg_value


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
            raise argparse.ArgumentTypeError("\n\n#ERROR : The parameter '--min-rdp-bootstrap' must be in format 'TAXONOMIC_LEVEL:MIN_BOOTSTRAP'.\n\n")
        output["rank"] = value.split(":")[0]
        output["value"] = value.split(":")[1]
        try:
            output["value"] = ratioParameter(output["value"])
        except:
            raise argparse.ArgumentTypeError("\n\n#ERROR : The value for the MIN_BOOTSTRAP in parameter '--min-rdp-bootstrap' must be between 0.0 and 1.0.\n\n")
        setattr(namespace, self.dest, output)


def impacted_obs_on_rdpBootstrap(input_biom, taxonomic_depth, min_bootstrap, impacted_file):
    """
    @summary: Writes the list of the observations with an insufficient bootstrap on the specified taxonomic rank.
    @param input_biom: [str] The path to the BIOM file to check.
    @param taxonomic_depth: [int] The taxonomic rank depth to check (example: 6 for Species in system "Domain, Phylum, Class, Order, Family, Genus, Species").
    @param min_bootstrap: [float] The observations with a value inferior to this threshold at the specified taxonomic depth are reported in the impacted file.
    @param impacted_file: [str] The path to the output file.
    @return n : number of filter OTU
    """
    n=0
    biom = BiomIO.from_json( input_biom )
    FH_impacted_file = open( impacted_file, "wt" )
    for observation in biom.get_observations():
        bootstrap = observation["metadata"]["rdp_bootstrap"]
        if issubclass(bootstrap.__class__, str):
            bootstrap = bootstrap.split(";")
        if bootstrap[taxonomic_depth] < min_bootstrap:
            FH_impacted_file.write( str(observation["id"]) + "\n" )
            n += 1
    FH_impacted_file.close()
    return n


def impacted_obs_on_blastMetrics( input_biom, tag, cmp_operator, threshold, impacted_file ):
    """
    @summary: Writes the list of the observations with no affiliations with sufficient blast value.
    @param input_biom: [str] The path to the BIOM file to check.
    @param tag: [str] The metadata checked.
    @param cmp_operator: [str] The operator use in comparison (tag_value ">=" thresold or tag_value "<=" thresold ).
    @param threshold: [float] The limit for the tag value.
    @param impacted_file: [str] The path to the output file.
    @return n : number of filtered OTU
    """
    n = 0
    valid_operators = {
        ">=": operator.__ge__,
        "<=": operator.__le__
    }
    cmp_func = valid_operators[cmp_operator]
    biom = BiomIO.from_json( input_biom )
    FH_impacted_file = open( impacted_file, "wt" )
    for observation in biom.get_observations():
        alignments = observation["metadata"]["blast_affiliations"]
        is_discarded = True
        for current_alignment in alignments:
            if cmp_func(float(current_alignment[tag]), threshold):
                is_discarded = False
        if is_discarded:
            FH_impacted_file.write( str(observation["id"]) + "\n" )
            n += 1

    FH_impacted_file.close()
    return n

def get_tax_consensus( taxonomies ):
    """
    @summary: Returns a consensus taxonomy from list of taxonomies.
    @param taxonomies: [list] The taxonomies to process. Each taxonomy is a list of rank taxon.
    @return: [list] The consensus taxonomy. The ambiguous ranks are replaced by "Multi-affiliation".
    @note:
        taxonomies = [ ["Bacteria", "Proteobacteria", "Gamma Proteobacteria", "Enterobacteriales"],
                       ["Bacteria", "Proteobacteria", "Beta Proteobacteria", "Methylophilales"] ]
        return = ["Bacteria", "Proteobacteria", "Multi-affiliation", "Multi-affiliation"]
    """
    consensus = list()
    if len(taxonomies) != 0:
        consensus = copy.copy(taxonomies[0])
    for curr_taxonomy in taxonomies[1:]:
        for rank, taxon in enumerate(curr_taxonomy):
            if consensus[rank] != "Multi-affiliation" and consensus[rank] != taxon:
                consensus[rank] = "Multi-affiliation"
    # Clean case with same taxon name in different branches:
    #      with taxonomies = [["A", "B", "C"], ["A", "L", "C"]]
    #      consensus is ["A", "Multi-affiliation", "C"] but must be ["A", "Multi-affiliation", "Multi-affiliation"]
    ancestor = ""
    for rank in range(len(consensus)):
        if ancestor == "Multi-affiliation":
            consensus[rank] = "Multi-affiliation"
        ancestor = consensus[rank]
    return consensus

def impacted_obs_by_undesired_taxon(input_biom, undesired_taxon_list, in_all_or_in_consensus, biom_out, impacted_file):
    """
    @summary : write the list of observation with affiliations including undesired taxon.
    @param input_biom: [str] The path to the BIOM file to check.
    @param undesired_taxon_list: [list] list of string to look for
    @param in_all_or_in_consensus: [bool] if True, one taxon_ignored must be in the consensus or all affiliation must contain one of the taxon ignored
    @param biom_out: [str] path to biom with removed undesired taxonomy
    @param impacted_file: [str] The path to the output file.
    @return n : number of filtered OTU
    """
    n = 0
    biom = BiomIO.from_json( input_biom )
    FH_impacted_file = open( impacted_file, "wt" )

    for observation in biom.get_observations():
        if not "comment" in observation['metadata'] :
            observation['metadata']['comment'] = list()

        # update blast_affiliations without ignored taxon and recompute de blast_taxonomy
        new_blast_affi = list()
        for affiliation in observation['metadata']['blast_affiliations']:
            keep=True
            for t in undesired_taxon_list:
                regexp = re.compile(t)
                if regexp.search(";".join(affiliation["taxonomy"])):
                    keep=False
            #~ if not any(t in ";".join(affiliation["taxonomy"]) for t in undesired_taxon_list):
            if keep:
                new_blast_affi.append(affiliation)

        # if some affi are masked, update blast_affiliations and blast_taxonomy
        if len(new_blast_affi) != len(observation['metadata']['blast_affiliations']):
            observation['metadata']['blast_affiliations'] = new_blast_affi
            new_consensus = get_tax_consensus([affi['taxonomy'] for affi in new_blast_affi] )
            
            # delete mode if all affiliations belong to one of undesired taxon
            if in_all_or_in_consensus and len(new_blast_affi) == 0 : 
                n += 1
                FH_impacted_file.write( str(observation["id"]) + "\n" )
            
            # masking mode if the new consensus is changed because of ignoring undesired taxon
            elif not in_all_or_in_consensus and new_consensus != observation['metadata']['blast_taxonomy']:
                n += 1
                FH_impacted_file.write( str(observation["id"]) + "\n" )

            # add comment
            if len(new_consensus) == 0:
                observation['metadata']['comment'].append("blast_tax_masked")
            elif not new_consensus == observation['metadata']['blast_taxonomy']:
                observation['metadata']['comment'].append("blast_tax_changed")

            observation['metadata']['blast_taxonomy'] = new_consensus

    BiomIO.write( biom_out, biom )
    return n


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

def list_from_file(in_file):
    uniq = list()
    FH_current_file = open( in_file, 'rt' )
    for line in FH_current_file:
        if not line.strip() in uniq:
            uniq.append(line.strip())
    FH_current_file.close()

    return uniq

def uniq_from_files_lists( in_files ):
    """
    @summary: Returns an list without duplicated elements from several list files.
    @param in_files: [list] The list of files paths. Each file contains a list.
    @return: [list] The list without duplicated elements.
    """
    uniq = list()
    for current_file in in_files:
        l = list_from_file(current_file)
        uniq = list(set(uniq + l ))

    return uniq

def mask_observation(rdp_clusters_discards, blast_clusters_discards, input_biom, output_biom):
    """
    @summary : mask either rdp affiliations and/or blast affiliations
    @param rdp_clusters_discards : [list] of clusters whith rdp affiliations to mask
    @param blast_clusters_discards : [list] of clusters whith blast consensus affiliations to mask
    @param input_biom : [str] Path to input biom file
    @param input_biom : [str] Path to output biom file with affiliations masked
    """

    biom = BiomIO.from_json(input_biom)
    for observation in biom.get_observations():
        if 'comment' not in observation['metadata'] :
            observation['metadata']['comment'] = list()
        # remove rdp taxonomic metadata
        if rdp_clusters_discards is not None and observation['id'] in rdp_clusters_discards:
            observation['metadata']['comment'].append("rdp_tax_masked")

            if issubclass( observation['metadata']["rdp_taxonomy"].__class__, str):
                observation['metadata']["rdp_taxonomy"] = ""
                observation['metadata']["rdp_bootstrap"] = ""
            elif issubclass( observation['metadata']["rdp_taxonomy"].__class__, list):
                observation['metadata']["rdp_taxonomy"] = list()
                observation['metadata']["rdp_bootstrap"] = list()

        # remove blast metadata
        if observation['id'] in blast_clusters_discards:
            observation['metadata']['comment'].append("blast_tax_masked")
            observation['metadata']["blast_affiliations"] = list()
            observation['metadata']["blast_taxonomy"] = list()

    BiomIO.write( output_biom, biom )

def write_impact( params, discards, tmp_biom, impacted_file, impacted_multihit ):
    """
    @summary: Writes the list of observations removed by each filter.
    @param params: [Namespace] to recover input bam file, and filter threshold
    @param discards: [dict] By filter the path of the file that contains the list of the removed observations.
    @param impacted_file: [str] The path to the abundance tsv file with filter annotation.
    @param impacted_multihit : [str] The path to the multihit tsv file for impacted OTU.
    """
    dict_OTU_discards = dict()

    # store reason why OTU are impacted
    for filter in discards:
        label=""
        if "min_rdp_bootstrap" in discards[filter]:
            label = params.min_rdp_bootstrap["rank"] + "_rdp_boostrap_lt_" + str(params.min_rdp_bootstrap["value"])
        elif "min_blast_length" in discards[filter]:
            label = "blast_len_lt_" + str(params.min_blast_length)
        elif "max_blast_evalue" in discards[filter]:
            label = "blast_evalue_gt_" + str(params.max_blast_evalue )
        elif "min_blast_identity" in discards[filter]:
            label = "blast_id_lt_" + str(params.min_blast_identity )
        elif "min_blast_coverage" in discards[filter]:
            label = "blast_cov_lt_" + str(params.min_blast_coverage )
        elif "taxon_ignored" in discards[filter]:
            label = "undesired_tax_in_blast"

        FH_disc = open(discards[filter],'rt')

        for line in FH_disc:
            otu_name = line.strip()
            if not otu_name in dict_OTU_discards:
                dict_OTU_discards[otu_name] = label
            else:
                dict_OTU_discards[otu_name] += ";"+label

    # Create a biom abundance file with impacted OTU only
    biom_In = BiomIO.from_json(params.input_biom)
    biom_Out = Biom( generated_by='FROGS_affiliation_filters', matrix_type="sparse" )

    # Add samples
    for sample_name in biom_In.get_samples_names():
        biom_Out.add_sample( sample_name )

    # Add impacted OTU
    for observation in biom_In.get_observations():
        if observation['id'] in dict_OTU_discards:
            if not 'comment' in observation['metadata']:
                observation['metadata']['comment'] = ""
            observation['metadata']['comment'] = dict_OTU_discards[observation['id']]
            if not 'affi_filter_statut' in observation['metadata']:
                observation['metadata']['affi_filter_statut'] = ""
            observation['metadata']['affi_filter_statut'] = 'affiliation_masked' if params.mask else 'OTU_deleted'
            biom_Out.add_observation( observation['id'], observation['metadata'] )
            # Add count
            for sample_name in biom_In.get_samples_names():
                if biom_In.get_count(observation['id'],sample_name) > 0:
                    biom_Out.add_count(observation['id'],sample_name, biom_In.get_count(observation['id'],sample_name))

    # Write
    BiomIO.write( tmp_biom, biom_Out )

    # Convert Biom in TSV
    Biom_to_tsv(tmp_biom, impacted_file, impacted_multihit).submit(params.log_file)




def write_summary( summary_file, input_biom, output_biom, discards, params ):
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

    # Number of taxon in each rank before filters
    in_biom = BiomIO.from_json( input_biom )
    taxon_lost = {'Blast': { rank : list() for rank in range(len(params.taxonomic_ranks))}}
    if in_biom.has_observation_metadata('rdp_taxonomy'):
        taxon_lost['RDP'] = {rank : list() for rank in range(len(params.taxonomic_ranks))}
    # for blast affiliation, consensus maybe updated, so lost of multiaffiliation should not be considered as taxon lost
    taxon_common = { rank : list() for rank in range(len(params.taxonomic_ranks))}
    taxon_gain = { rank : list() for rank in range(len(params.taxonomic_ranks))}

    for observation_name in in_biom.get_observations_names():
        global_results['nb_clstr_ini'] += 1
        global_results['nb_seq_ini'] += in_biom.get_observation_count( observation_name )
        # track rdp taxon
        if 'RDP' in taxon_lost:
            rdp_taxonomy = in_biom.get_observation_metadata(observation_name)['rdp_taxonomy']
            if issubclass(rdp_taxonomy.__class__,str):
                rdp_taxonomy = rdp_taxonomy.split(';')
            for i in range(len(params.taxonomic_ranks)):
                if ';'.join(rdp_taxonomy[:i+1]) not in taxon_lost['RDP'][i]:
                    taxon_lost['RDP'][i].append(';'.join(rdp_taxonomy[:i+1]))
        # track blast taxon
        blast_taxonomy = in_biom.get_observation_metadata(observation_name)['blast_taxonomy']
        if issubclass(blast_taxonomy.__class__,str):
            blast_taxonomy = blast_taxonomy.split(';')
        for i in range(len(params.taxonomic_ranks)):
            if ';'.join(blast_taxonomy[:i+1]) not in taxon_lost['Blast'][i]:
                taxon_lost['Blast'][i].append(';'.join(blast_taxonomy[:i+1]))

    for sample_name in in_biom.get_samples_names():
        samples_results[sample_name] = {
            'initial': sum( 1 for x in in_biom.get_observations_by_sample(sample_name) ),
            'filtered': dict(),
            'kept': 0
        }


    # By sample and by filters
    filters_intersections = dict()
    for filter in list(discards.keys()):
        FH_filter = open( discards[filter] )
        for line in FH_filter:
            observation_name = line.strip()
            if observation_name not in filters_intersections:
                filters_intersections[observation_name] = dict()
            filters_intersections[observation_name][filter] = 1
        FH_filter.close()
    for observation_name in list(filters_intersections.keys()):
        # Removed intersection
        intersections_key = "--@@--".join(sorted( filters_intersections[observation_name].keys() ))
        if intersections_key not in filters_results:
            filters_results[intersections_key] = {
                'filters': list(filters_intersections[observation_name].keys()),
                'count': 0
            }
        filters_results[intersections_key]['count'] += 1

        # Filters by samples
        for sample in in_biom.get_samples_by_observation(observation_name):
            for filter in filters_intersections[observation_name]:
                if filter not in samples_results[sample['id']]['filtered']:
                    samples_results[sample['id']]['filtered'][filter] = 0
                samples_results[sample['id']]['filtered'][filter] += 1
    del in_biom

    # Global after filters
    out_biom = BiomIO.from_json( output_biom )
    for observation_name in out_biom.get_observations_names():
        global_results['nb_clstr_kept'] += 1
        global_results['nb_seq_kept'] += out_biom.get_observation_count( observation_name )

        # track lost rdp taxon
        if 'RDP' in taxon_lost:
            rdp_taxonomy = out_biom.get_observation_metadata(observation_name)['rdp_taxonomy']
            if issubclass(rdp_taxonomy.__class__,str):
                rdp_taxonomy = rdp_taxonomy.split(';')
            for i in range(len(params.taxonomic_ranks)):
                if ';'.join(rdp_taxonomy[:i+1]) in taxon_lost['RDP'][i]:
                    taxon_lost['RDP'][i].remove(';'.join(rdp_taxonomy[:i+1]))

        # track blast taxon
        blast_taxonomy = out_biom.get_observation_metadata(observation_name)['blast_taxonomy']
        if issubclass(blast_taxonomy.__class__,str):
            blast_taxonomy = blast_taxonomy.split(';')
        for i in range(len(params.taxonomic_ranks)):
            if ';'.join(blast_taxonomy[:i+1]) in taxon_lost['Blast'][i]:
                taxon_lost['Blast'][i].remove(';'.join(blast_taxonomy[:i+1]))
                taxon_common[i].append(';'.join(blast_taxonomy[:i+1]))
            elif ';'.join(blast_taxonomy[:i+1]) not in taxon_lost['Blast'][i] + taxon_common[i] + taxon_gain[i]:
                taxon_gain[i].append(';'.join(blast_taxonomy[:i+1]))

    for sample_name in out_biom.get_samples_names():
        samples_results[sample_name]['kept'] = sum( 1 for x in out_biom.get_observations_by_sample(sample_name) )
    del out_biom

    #print(taxon_lost['Blast'][6])
    #print(taxon_gain[6])
    # store number instead of list of taxon lost
    for i in range(len(params.taxonomic_ranks)):
        taxon_lost['Blast'][params.taxonomic_ranks[i]] = len(taxon_lost['Blast'][i]) - len(taxon_gain[i])
        taxon_lost['Blast'].pop(i)
        if 'RDP' in taxon_lost:
            taxon_lost['RDP'][params.taxonomic_ranks[i]] = len(taxon_lost['RDP'][i])
            taxon_lost['RDP'].pop(i)

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "affiliation_filters_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###PORCESSED_FILTERS###" in line:
            line = line.replace( "###PORCESSED_FILTERS###", json.dumps([filter for filter in discards]) )
        elif "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###GLOBAL_TAXON_LIST###" in line:
            line = line.replace( "###GLOBAL_TAXON_LIST###", json.dumps(params.taxonomic_ranks) )
        elif "###GLOBAL_TAXON_LOST###" in line:
            line = line.replace( "###GLOBAL_TAXON_LOST###", json.dumps(taxon_lost) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
        elif "###FILTERS_RESULTS###" in line:
            line = line.replace( "###FILTERS_RESULTS###", json.dumps(list(filters_results.values())) )
        elif "Draw a Venn to see which OTUs had been deleted by the filters chosen (Maximum 6 options): " in line and params.mask:
            line = "Draw a Venn to see which OTUs had its taxonomy masked by the filters chosen (Maximum 6 options): "
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()


def process( args ):
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )

    biom_in = args.input_biom
    try:
        discards = dict() # by filter the discard file path
        Logger.static_write(args.log_file, "Identify OTU with :\n")
        if args.min_rdp_bootstrap is not None:
            label = "RDP bootstrap for " + args.min_rdp_bootstrap["rank"] + " < " + str(args.min_rdp_bootstrap["value"])
            discards[label] = tmpFiles.add( "min_rdp_bootstrap" )
            n = impacted_obs_on_rdpBootstrap( biom_in, args.rdp_taxonomy_ranks.index(args.min_rdp_bootstrap["rank"]), args.min_rdp_bootstrap["value"], discards[label] )
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")

        if args.min_blast_length is not None:
            label = "All blast length < " + str(args.min_blast_length)
            discards[label] = tmpFiles.add( "min_blast_length" )
            n = impacted_obs_on_blastMetrics( biom_in, "aln_length", ">=", args.min_blast_length, discards[label] )
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")

        if args.max_blast_evalue is not None:
            label = "All blast evalue > " + str(args.max_blast_evalue)
            discards[label] = tmpFiles.add( "max_blast_evalue" )
            n = impacted_obs_on_blastMetrics( biom_in, "evalue", "<=", args.max_blast_evalue, discards[label] )
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")

        if args.min_blast_identity is not None:
            label = "All blast identity < " + str(args.min_blast_identity)
            discards[label] = tmpFiles.add( "min_blast_identity" )
            n = impacted_obs_on_blastMetrics( biom_in, "perc_identity", ">=", 100*args.min_blast_identity, discards[label] )
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")

        if args.min_blast_coverage is not None:
            label = "All blast coverage < " + str(args.min_blast_coverage)
            discards[label] = tmpFiles.add( "min_blast_coverage" )
            n = impacted_obs_on_blastMetrics( biom_in, "perc_query_coverage", ">=", 100*args.min_blast_coverage, discards[label] )
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")
        
        if args.taxon_ignored is not None:
            label = "All blast taxonomies or consensus taxonomy belong to undesired taxon: " + " / ".join(args.taxon_ignored)
            if args.mask:
                label = "Some blast taxonomies belong to undesired taxon: " + " / ".join(args.taxon_ignored)
            discards[label] = tmpFiles.add( "taxon_ignored" )
            biom_taxIgnored = tmpFiles.add( "taxon_ignored.biom" )
            n = 0
            if args.delete:
                n = impacted_obs_by_undesired_taxon(biom_in, args.taxon_ignored, True, biom_taxIgnored, discards[label])
            else:
                n = impacted_obs_by_undesired_taxon(biom_in, args.taxon_ignored, False, biom_taxIgnored, discards[label])
            Logger.static_write(args.log_file, "\t- " + label + " : "+ str(n) + "\n")
            biom_in = biom_taxIgnored


        if args.delete:
            clusters_discarded = uniq_from_files_lists( [discards[filter] for filter in discards] )
            remove_observations( clusters_discarded, biom_in, args.output_biom )
            update_fasta_log = tmpFiles.add( "update_fasta_log.txt" )
            UpdateFasta( args.output_biom, args.input_fasta, args.output_fasta, update_fasta_log ).submit( args.log_file )
        
        if args.mask:
            rdp_clusters_discards = list()
            if args.min_rdp_bootstrap is not None and "RDP bootstrap for " + args.min_rdp_bootstrap["rank"] + " < " + str(args.min_rdp_bootstrap["value"]) in discards: 
                rdp_clusters_discards = list_from_file(discards["RDP bootstrap for " + args.min_rdp_bootstrap["rank"] + " < " + str(args.min_rdp_bootstrap["value"])])
            blast_clusters_discards = uniq_from_files_lists( [discards[filter] for filter in discards if "All blast" in filter and not "undesired taxon" in filter] )
            mask_observation(rdp_clusters_discards, blast_clusters_discards, biom_in, args.output_biom)

        # Writes outputs
        impacted_biom = tmpFiles.add("impacted.biom")
        write_impact( args, discards, impacted_biom, args.impacted, args.impacted_multihit)
        write_summary( args.summary, args.input_biom, args.output_biom, discards, args )

    finally:
        if not args.debug : 
            tmpFiles.deleteAll()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == '__main__':
    # Parameters
    parser = argparse.ArgumentParser(description='Filters an abundance BIOM file on affiliations metrics')
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '--taxonomic-ranks', nargs='+', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels used in the metadata taxonomy. [Default: %(default)s]' )
    #     Filters behavior
    group_filter_bh = parser.add_argument_group( 'Filters behavior' )
    group_exclusion_filter_bh = group_filter_bh.add_mutually_exclusive_group()
    group_exclusion_filter_bh.add_argument('-m','--mask', default=False, action='store_true', help="If affiliations do not respect one of the filter they are replaced by NA (mutually exclusive with --delete)")
    group_exclusion_filter_bh.add_argument('-d','--delete', default=False, action='store_true', help="If affiliations do not respect one of the filter the entire OTU is deleted.(mutually exclusive with --mask)")
    #     Filters
    group_filter = parser.add_argument_group( 'Filters' )
    group_filter.add_argument( '--taxon-ignored', type=str, nargs='*', help="Taxon list to maks/delete")
    group_filter.add_argument( '-b', '--min-rdp-bootstrap', type=str, action=BootstrapParameter, metavar=("TAXONOMIC_LEVEL:MIN_BOOTSTRAP"), help="The minimal RDP bootstrap must be superior to this value (between 0 and 1)." )
    group_filter.add_argument( '-t', '--rdp-taxonomy-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels present in the reference databank. [Default: %(default)s]' )
    group_filter.add_argument( '-i', '--min-blast-identity', type=ratioParameter, help="The number corresponding to the blast percentage identity (between 0 and 1)." )
    group_filter.add_argument( '-c', '--min-blast-coverage', type=ratioParameter, help="The number corresponding to the blast percentage coverage (between 0 and 1)." )
    group_filter.add_argument( '-e', '--max-blast-evalue', type=float, help="The number corresponding to the blast e value (between 0 and 1).")
    group_filter.add_argument( '-l', '--min-blast-length', type=int, default=None, required=False, help="The number corresponding to the blast length." )
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help="The input BIOM file.")
    group_input.add_argument('--input-fasta', required=True, help="The input FASTA file.")
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default="affiliation_filters.biom", help="The BIOM file output. [Default: %(default)s]")
    group_output.add_argument('--output-fasta', default="affiliation_filters.fasta", help="The FASTA output file. [Default: %(default)s]")
    group_output.add_argument('--summary', default="affiliation_filters.html", help="The HTML file containing the graphs. [Default: %(default)s]")
    group_output.add_argument('--impacted', default="affiliation_filters_impacted_clusters.tsv", help="The abundance file that summarizes all the clusters impacted (deleted or with affiliations masked). [Default: %(default)s]")
    group_output.add_argument('--impacted-multihit', default="affiliation_filters_impacted_clusters_multihit.tsv", help="The multihit TSV file associated with impacted OTU. [Default: %(default)s]")
    group_output.add_argument('--log-file', default=sys.stdout, help='This output will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)


    # keep quote around each taxon ignored when printing command line into logfile
    cmd=list()
    taxon = False
    for arg in sys.argv:
        if arg == "--taxon-ignored" :
            taxon = True
            cmd.append(arg)
        elif taxon and arg.startswith("-"):
            taxon = False
        elif taxon:
            cmd.append('\"' + arg + '\"')
        if not taxon:
            cmd.append(arg)

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(cmd) + "\n\n")

    if not args.delete and not args.mask:
        raise argparse.ArgumentTypeError("\n\n#ERROR : You must precise if you want to mask affiliations of delete OTU with --mask or --delete options.\n\n")

    if args.min_rdp_bootstrap is None and args.min_blast_identity is None and args.min_blast_coverage is None and args.max_blast_evalue is None and args.min_blast_length is None and args.taxon_ignored is None:
        raise argparse.ArgumentTypeError( "\n\n#ERROR : At least one filter must be set to run " + os.path.basename(sys.argv[0]) + "\n\n")
    
    in_biom = BiomIO.from_json( args.input_biom )

    if args.min_rdp_bootstrap is not None:
        if not in_biom.has_observation_metadata("rdp_bootstrap"):
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'rdp_bootstrap'. You cannot use the parameter '--min-rdp-bootstrap' on this file.\n\n" )
        elif not args.min_rdp_bootstrap["rank"] in args.rdp_taxonomy_ranks:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The taxonomic rank choosen in '--min-rdp-bootstrap' must be in '--rdp-taxonomy-ranks' (" + ";".join(args.rdp_taxonomy_ranks) + ").\n\n" )
        else:
            nb_rank = 0
            for observation in in_biom.get_observations():
                if 'rdp_bootstrap' in observation['metadata'] and len(observation['metadata']['rdp_bootstrap']) > 0:
                    nb_rank = len(observation['metadata']['rdp_bootstrap'])
                    break
            if nb_rank != len(args.taxonomic_ranks):
                raise argparse.ArgumentTypeError("\n\n#ERROR : RDP taxonomic metadata is defined on " + str(nb_rank) + " and you precise " + str(len(args.taxonomic_ranks)) + " ranks name (see --taxonomic-ranks)\n\n")
    if not in_biom.has_observation_metadata("blast_affiliations"):
        if args.min_blast_length is not None:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-length' on this file.\n\n" )
        if args.max_blast_evalue is not None:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-evalue' on this file.\n\n" )
        if args.min_blast_identity is not None:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-identity' on this file.\n\n" )
        if args.min_blast_coverage is not None:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-coverage' on this file.\n\n" )
        if args.taxon_ignored is not None:
            raise argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--taxon-ignored' on this file.\n\n" )
    del in_biom

    # Process
    process( args )
