#!/usr/bin/env python3

__author__ = 'Katia Vidal - GENPHYSE & Frédéric Escudié - Genotoul/MIAT & Maria Bernard - SIGENAE/GABI'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
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
# CLASS
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
        out_files.append( tmp_files_manager.add(out_basename_pattern.replace('##RANK##', 'asv')) )

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
        
    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()        

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
        lower_quartile_idx = int(nb_samples/4)
        nb_seq = counts[lower_quartile_idx]
        # If lower quartile sample is empty
        if nb_seq == 0:
            idx = 1
            while (lower_quartile_idx + idx) < nb_samples and counts[lower_quartile_idx + idx] == 0:
                idx += 1
            if (lower_quartile_idx + idx) < nb_samples:
                nb_seq = counts[lower_quartile_idx + idx]
        step_size = int(nb_seq/nb_step)
        return max(1, step_size)


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
                      
    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()      
        
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
            raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : The parameter '--min-rdp-bootstrap' must be in format 'TAXONOMIC_LEVEL:MIN_BOOTSTRAP'.\n\n"))
        output["rank"] = value.split(":")[0]
        output["value"] = value.split(":")[1]
        try:
            output["value"] = ratioParameter(output["value"])
        except:
            raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : The value for the MIN_BOOTSTRAP in parameter '--min-rdp-bootstrap' must be between 0.0 and 1.0.\n\n"))
        setattr(namespace, self.dest, output)

class BIOM_to_TSV(Cmd):
    """
    @summary: Convert a biom file to a tabular file
    """
    def __init__(self, in_biom, in_fasta, out_tsv, out_multihit, out_log, header_only=False):
        """
        @param in_biom [str] : Path of the biom file to convert
        @param in_fasta [str] : Path of the fasta file associated with biom file
        @param out_tsv [str] : Path of the tabular file to create
        @param out_multihit [str] : Path of the detailed multiaffiliated ASV.
        @param out_log [str] : biom_to_tsv log file
        @param header_only [bool] : extract only header
        """
        opt = ''
        if header_only:
            opt = ' --header '
        Cmd.__init__( self,
                      'biom_to_tsv.py',
                      'Convert a biom file into tsv file',
                      " --input-biom " + in_biom + ' --input-fasta ' + in_fasta +' --output-tsv ' + out_tsv + ' --output-multi-affi ' + out_multihit + ' --log-file  ' + out_log + opt,
                      '--version' )

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').strip()

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
        return Cmd.get_version(self, 'stdout').strip()                      

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def ratioParameter( arg_value ):
    """
    @summary: Argparse type for ratio (float between 0 and 100).
    """
    float_arg_value = None
    try:
        float_arg_value = float(arg_value)
        if float_arg_value < 0.0 or float_arg_value > 100.0:
            raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : must be between 0 and 100.\n"))
    except:
        raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : must be between 0 and 100.\n"))
    return float_arg_value

def checkBlastTaxa( param ):
    if param is not None:
        temp = [ i.strip() for i in param if i.strip() != '']
        if len(temp) == 0:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : keep-blast-taxa or ignore-blast-taxa list need to be non empty strings.\n\n" ))
        else:
            if len(temp) != len(param):
                param = temp
                Logger.static_write(args.log_file, "WARNING : empty string in keep-blast-taxa or ignore-blast-taxa option have been removed, here is the updated list that will be take into account: \"" + '\" \"'.join(args.ignore_blast_taxa) + "\n\n")

def impacted_obs_on_rdpBootstrap(observation, taxonomic_depth, min_bootstrap):
    """
    @summary: check whether the observations is with an insufficient bootstrap on the specified taxonomic rank.
    @param observation: [obj] observation object.
    @param taxonomic_depth: [int] The taxonomic rank depth to check (example: 6 for Species in system "Domain, Phylum, Class, Order, Family, Genus, Species").
    @param min_bootstrap: [float] The bootstrap threshold.
    @return True or False
    """
    bootstrap = observation["metadata"]["rdp_bootstrap"]
    if issubclass(bootstrap.__class__, str):
        bootstrap = bootstrap.split(";")
    if bootstrap[taxonomic_depth] < min_bootstrap:
        return True
    else:
        return False

def impacted_blast_affi_on_blastMetrics( observation, tag, cmp_operator, threshold):
    """
    @summary: return blast_affiliation that respect the criteria
    @param observation: [object] observation object with a list of blast affiliations
    @param tag: [str] The blast affiliation metadata to check.
    @param cmp_operator: [str] The operator use in comparison (tag_value ">=" thresold or tag_value "<=" thresold ).
    @param threshold: [float] The limit for the tag value.
    @return blast affiliation filtered list
    """

    valid_operators = {
        ">=": operator.__ge__,
        "<=": operator.__le__
    }
    cmp_func = valid_operators[cmp_operator]
    blast_affiliations_out = dict()
    for idx,blast_affi in enumerate(observation['metadata']['blast_affiliations']):
        if blast_affi[tag] == "no data":
            sys.exit("\nError: " + observation["id"]+ " has a taxonomy but some of its metrics have 'no data'!\nHave you modified it manually? Please correct it.\n\n" + str(observation)+"\n\n")
        elif cmp_func(float(str(blast_affi[tag]).replace(',', '.')), threshold):
            blast_affiliations_out[idx] = blast_affi

    return blast_affiliations_out

def impacted_blast_affi_on_blastTaxonomy(observation, taxon_ignored, init_keep):
    """
    @summary: return blast affiliations whithout undesired taxon
    @param observation [obj] : observation object with list of blast affiliations
    @param taxon_ignored [list] : list of taxon to ignored (it may be partial terms)
    @param init_keep [boolean] : True if --ignore-blast-taxa, False if --keep-blast-taxa
    @return blast affiliations filtered list
    """
    blast_affiliations_out = dict()
    for idx,blast_affi in enumerate(observation['metadata']['blast_affiliations'] ):
        # True if ignore_blast_taxa, False if keep_blast_taxa
        keep = init_keep
        for t in taxon_ignored:
            regexp = re.compile(t)
            if regexp.search(";".join(blast_affi["taxonomy"])):
                if keep:
                    keep=False
                else:
                    keep=True
        if keep:
            blast_affiliations_out[idx] = blast_affi

    return blast_affiliations_out

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

    if len(consensus) == 0 :
        consensus = None
    return consensus



def update_blast_metadata(metadata_dict, kept_affiliaitons):
    """
    @summary : return impact status
    @param metadata_dict [dict] : original observation metadata dictionnary
    @param kept_affiliations [dict] : indexed blast affiliations to keep by filter
    @return a filtering status True (if all affiliation are removed), False otherwise
    """
    
    # compare the number of filtering criteria to the number of time a taxonomy is kept.
    nb_criteria = len(kept_affiliaitons)
    nb_index = dict()
    for criteria in kept_affiliaitons:
        for i in kept_affiliaitons[criteria]:
            if i not in nb_index:
                nb_index[i] = 0
            nb_index[i] += 1

    keys = list(nb_index.keys())
    keys.sort(reverse=True)
    
    # keep affiliation present in all filtered list
    for index in range(len(metadata_dict['blast_affiliations'])-1 , -1, -1) :
        if index not in nb_index or nb_index[index] < nb_criteria:
            metadata_dict['blast_affiliations'].pop(index)

    # update consensus taxonomy
    metadata_dict['blast_taxonomy'] = get_tax_consensus([affi['taxonomy'] for affi in metadata_dict['blast_affiliations']] )

    # return impacting status
    if metadata_dict['blast_taxonomy'] is None :
        metadata_dict['blast_affiliations'] = []
        return True
    else:
        return False


def add_obs(input_biom, observation_name, observation_metadata, output_biom):
    """
    @summary : add observation and sample count from an input biom to an output biom with updated metadata
    @param input_biom [Obj] : biom object corresponding to original biom file
    @parma observation_name [str] : name of the observation to add
    @param observation_metadata [dict] : dictionnary of observation metadata to add
    @param output biom [Obj] : biom object corresponding to output biom file
    """
    output_biom.add_observation( observation_name, observation_metadata )
    # Add count
    for sample_name in input_biom.get_samples_names():
        if input_biom.get_count(observation_name,sample_name) > 0:
            output_biom.add_count(observation_name,sample_name, input_biom.get_count(observation_name,sample_name))


def get_uniq_tax ( indexed_affiliations):
    """
    @summary : return uniq taxonomy list
    @param indexed_affiliations [dict] : dictionnary of indexed blast_affiliations {0:{"perc_identity": 100.0, "taxonomy": ["Bacteria", "Bacteroidetes", "Flavobacteriia", "Flavobacteriales", "Flavobacteriaceae", "Pibocella", "Pibocella ponti"], "evalue": "0.0", "aln_length": 421, "perc_query_coverage": 100.0, "subject": "AY576654.1.1447"}}
    """
    uniq_tax=list()
    for idx in indexed_affiliations :
        tax = indexed_affiliations[idx]['taxonomy']
        if issubclass( tax.__class__, list):
            tax = ";".join(tax)
        if not tax in uniq_tax:
            uniq_tax.append(tax)
    return uniq_tax

def filter_biom(in_biom_file, impacted_file, output_file, params):
    """
    @summary : parse in_biom and delete ASV or mask affiliation whether they do not respect filter expressed in params
    @param in_biom [str] : Path to input biom file
    @param impacted_file [str] : Path to impacted biom file to store all impacted (deleted / masked or updated) ASV
    @param output_file [str] : Path to output biom file
    @param params [NameSpace] : taxonomical filtering criteria
    @return impacted_dict to summarize the results
    """
    # init
    impacted_dict=dict()
    in_biom = BiomIO.from_json( in_biom_file )

    # create output biom.
    impacted_biom = Biom( generated_by='FROGS_affiliation_filters', matrix_type="sparse" )
    out_biom = Biom( generated_by='FROGS_affiliation_filters', matrix_type="sparse" )
    # Add samples
    for sample_name in in_biom.get_samples_names():
        impacted_biom.add_sample( sample_name )
        out_biom.add_sample( sample_name )

    # parse input
    for observation in in_biom.get_observations():
        if not 'comment' in observation['metadata']:
            observation['metadata']['comment']=list()
        # in old version of FROGS comment may be a string
        elif issubclass(observation['metadata']['comment'].__class__, str):
            observation['metadata']['comment'] = observation['metadata']['comment'].split(";")

        # filter on RDP boostrap criteria
        filter_on_rdpBootstrap = True
        if params.min_rdp_bootstrap is not None and observation['metadata']['rdp_taxonomy'] is not None:

            # add rdp bootstrap criteria
            label = "RDP bootstrap for " + params.min_rdp_bootstrap["rank"] + " < " + str(params.min_rdp_bootstrap["value"])
            
            filter_on_rdpBootstrap = impacted_obs_on_rdpBootstrap(observation, params.taxonomic_ranks.index(params.min_rdp_bootstrap["rank"]), params.min_rdp_bootstrap["value"])
            if filter_on_rdpBootstrap :
                if not label in impacted_dict:
                    impacted_dict[label] = list()
                impacted_dict[label].append(observation['id'])
                observation['metadata']['comment'].append(params.min_rdp_bootstrap["rank"] + "_rdp_boostrap_lt_" + str(params.min_rdp_bootstrap["value"]))
        else : 
            filter_on_rdpBootstrap = False

        # to store by criteria valid blast affiliations 
        filter_on_blastCriteria = False
        out_blast_affiliations = dict()
        # to check if criteria has an impact
        tax_in=list()
        if observation['metadata']['blast_affiliations'] :
            for blast_affi in observation['metadata']['blast_affiliations']:
                tax = blast_affi['taxonomy']
                if issubclass( blast_affi["taxonomy"].__class__, list):
                    tax = ";".join(blast_affi["taxonomy"])
                if not tax in tax_in:
                    tax_in.append(tax)

            # add blast length criteria
            if params.min_blast_length:
                label = "Blast length < " + str(params.min_blast_length)
                out_blast_affiliations['filter_on_len'] = impacted_blast_affi_on_blastMetrics(observation, "aln_length", ">=", params.min_blast_length)
                uniq_tax = get_uniq_tax(out_blast_affiliations['filter_on_len'])
                if len(uniq_tax) != len(tax_in):
                    observation['metadata']['comment'].append("blast_len_lt_" + str(params.min_blast_length))
                    if not label in impacted_dict:
                        impacted_dict[label] = list()
                    impacted_dict[label].append(observation['id'])
                elif len(uniq_tax) == len(tax_in):
                    out_blast_affiliations.pop('filter_on_len')  

            # add blast evalue criteria
            if params.max_blast_evalue:
                label = "Blast evalue > " + str(args.max_blast_evalue)
                out_blast_affiliations['filter_on_evalue'] = impacted_blast_affi_on_blastMetrics(observation, "evalue", "<=", params.max_blast_evalue)
                uniq_tax = get_uniq_tax(out_blast_affiliations['filter_on_evalue'])
                if len(uniq_tax) != len(tax_in):
                    observation['metadata']['comment'].append("blast_eval_gt_" + str(params.max_blast_evalue))
                    if not label in impacted_dict:
                        impacted_dict[label] = list()
                    impacted_dict[label].append(observation['id'])
                elif len(uniq_tax) == len(tax_in):
                    out_blast_affiliations.pop('filter_on_evalue')  

            # add blast identity criteria
            if params.min_blast_identity:
                label = "Blast identity < " + str(args.min_blast_identity)
                out_blast_affiliations['filter_on_identity'] = impacted_blast_affi_on_blastMetrics(observation, "perc_identity", ">=", params.min_blast_identity)
                uniq_tax = get_uniq_tax(out_blast_affiliations['filter_on_identity'])
                if len(uniq_tax) != len(tax_in):
                    observation['metadata']['comment'].append("blast_identity_lt_" + str(params.min_blast_identity))
                    if not label in impacted_dict:
                        impacted_dict[label] = list()
                    impacted_dict[label].append(observation['id'])
                elif len(uniq_tax) == len(tax_in):
                    out_blast_affiliations.pop('filter_on_identity')         

            # add blast coverage criteria
            if params.min_blast_coverage:
                label = "Blast coverage < " + str(args.min_blast_coverage)
                out_blast_affiliations['filter_on_coverage'] = impacted_blast_affi_on_blastMetrics(observation, "perc_query_coverage", ">=", params.min_blast_coverage)
                uniq_tax = get_uniq_tax(out_blast_affiliations['filter_on_coverage'])
                if len(uniq_tax) != len(tax_in):
                    observation['metadata']['comment'].append("blast_coverage_lt_" + str(params.min_blast_coverage))
                    if not label in impacted_dict:
                        impacted_dict[label] = list()
                    impacted_dict[label].append(observation['id'])
                elif len(uniq_tax) == len(tax_in):
                    out_blast_affiliations.pop('filter_on_coverage')

            # add blast taxon to ignore criteria
            if params.ignore_blast_taxa or params.keep_blast_taxa:
                if params.ignore_blast_taxa:
                    param = params.ignore_blast_taxa
                    keep = True
                    label = "Blast taxonomies belong to undesired taxon: " + " / ".join(param)
                else:
                    param = params.keep_blast_taxa
                    keep = False
                    label = "Blast taxonomies belong to desired taxon: " + " / ".join(param)
                out_blast_affiliations['filter_on_taxonIgnored'] = impacted_blast_affi_on_blastTaxonomy(observation, param, keep)
                uniq_tax = get_uniq_tax(out_blast_affiliations['filter_on_taxonIgnored'])
                if len(uniq_tax) != len(tax_in):
                    observation['metadata']['comment'].append("undesired_tax_in_blast")
                    if not label in impacted_dict:
                        impacted_dict[label] = list()
                    impacted_dict[label].append(observation['id'])
                elif len(uniq_tax) == len(tax_in):
                    out_blast_affiliations.pop('filter_on_taxonIgnored')

        elif params.min_blast_length or params.max_blast_evalue or params.min_blast_identity or params.min_blast_coverage or params.ignore_blast_taxa or params.keep_blast_taxa and not observation['metadata']['blast_affiliations']:
            label = "Blast missing affiliations"
            if not label in impacted_dict:
                impacted_dict[label] = list()
            impacted_dict[label].append(observation['id'])
            filter_on_blastCriteria = True
            
        # update blast_affiliation dictionnary and compute blast filtering status : True if all affiliations are removed else False
        metadata_out = copy.deepcopy(observation['metadata'])
        if len(out_blast_affiliations) > 0:
            filter_on_blastCriteria = update_blast_metadata(metadata_out, out_blast_affiliations)

        # keep impacting criteria only if filter_on_blastCriteria is True
        if not filter_on_blastCriteria : 
            for impact in impacted_dict:
                if impact.startswith('Blast') and observation['id'] in impacted_dict[impact]:
                    impacted_dict[impact].remove(observation['id'])

        # write observation in impacted biom as the orignal but with additionnal status metadata corresponding to 
        # the type of filtering (ASV_deleted/Affiliation_masked/Blast_taxonomy_changed) and/or in output biom file whithout affiliation that do not respect one of the criteria
        if params.delete :
            # set status to ASV_deleted, blast_taxonomy changed or None
            if filter_on_rdpBootstrap or filter_on_blastCriteria:
                observation['metadata']['status'] = 'ASV_deleted'
                add_obs(in_biom, observation['id'], observation['metadata'], impacted_biom)
                
            # if change in blast consensus taxonomy, write ASV in both output biom and impacted biom
            else:
                add_obs(in_biom, observation['id'], metadata_out, out_biom)
                # if blast consensus taxonomy changed, store also in impacted biom the original ASV.
                if observation['metadata']['blast_taxonomy'] != metadata_out['blast_taxonomy']:
                    if not 'Blast_taxonomy_changed' in impacted_dict:
                        impacted_dict['Blast_taxonomy_changed'] = list()
                    impacted_dict["Blast_taxonomy_changed"].append(observation['id'])
                    observation['metadata']['status'] = "Blast_taxonomy_changed"
                    add_obs(in_biom, observation['id'], observation['metadata'], impacted_biom)
            
        elif params.mask :

            status = None
            # mask RDP taxonomy
            if filter_on_rdpBootstrap:
                metadata_out['rdp_taxonomy'] =  None
                metadata_out['rdp_bootstrap'] =  None

            # compute status either None, or Affiliation_masked +/- Blast_taxonomy_changed
            if filter_on_rdpBootstrap or filter_on_blastCriteria:
                status = "Affiliation_masked"
            if metadata_out['blast_taxonomy'] and observation['metadata']['blast_taxonomy'] != metadata_out['blast_taxonomy']:
                if status :
                    status += ";Blast_taxonomy_changed"
                else:
                    status = "Blast_taxonomy_changed"

            # write ASV with possible metadata updated into output biom
            add_obs(in_biom, observation['id'], metadata_out, out_biom)
            # add status metadata and write ASV in impacted biom if any change
            if status :
                observation['metadata']['status'] = status
                add_obs(in_biom, observation['id'], observation['metadata'], impacted_biom)
                # add ASV name in list of ASV impacted, ASV have not been mask but only updated
                if status == "Blast_taxonomy_changed":
                    if status not in impacted_dict:
                        impacted_dict[status] = list()
                    impacted_dict[status].append(observation['id'])

    # write output and impacted biom
    BiomIO.write( output_file, out_biom )
    BiomIO.write( impacted_file, impacted_biom )

    return impacted_dict

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
        identity = 0
        coverage = 0
        if multiple_tag is not None:
            if multiple_tag in observation_metadata and observation_metadata[multiple_tag] is not None and len(observation_metadata[multiple_tag]) > 0:
                identity = observation_metadata[multiple_tag][0][identity_tag]
                coverage = observation_metadata[multiple_tag][0][coverage_tag]
        else:
            if identity_tag in observation_metadata and coverage_tag in observation_metadata:
                identity = observation_metadata[identity_tag]
                coverage = observation_metadata[coverage_tag]
        if identity not in aln_results_hash:
            aln_results_hash[identity] = dict()
        if coverage not in aln_results_hash[identity]:
            aln_results_hash[identity][coverage] = {
                "clstr": 0,
                "seq": 0
            }
        aln_results_hash[identity][coverage]["clstr"] += 1
        aln_results_hash[identity][coverage]["seq"] += biom.get_observation_count( observation['id'] )
    
    for ident in list(aln_results_hash.keys()):
        for cover in list(aln_results_hash[ident].keys()):
            aln_results.append([
                ident,
                cover,
                aln_results_hash[ident][cover]
            ])
    del biom
    return aln_results

def write_summary( summary_file, input_biom, output_biom, discards, tree_count_file, tree_ids_file, rarefaction_files, rarefaction_ranks, params ):
    """
    @summary: Writes the process summary.
    @param summary_file: [str] The path to the output HTML file.
    @param input_biom: [str] The path to the BIOM before program execution.
    @param output_biom: [str] The path to the BIOM after program execution.
    @param discards: [dict] By filter the list of the deleted/masked observations.
    """
    
    # Get rarefaction data
    rarefaction_step_size = None
    rarefaction = None
    biom = BiomIO.from_json( input_biom )
    rarefaction_ranks.append('ASVs')
    for rank_idx, current_file in enumerate(rarefaction_files):
        rank = rarefaction_ranks[rank_idx]
        FH_rarefaction = open( current_file )
        for line in FH_rarefaction:
            fields = list(map(str.strip, line.split("\t")))
            if line.startswith('#'):
                samples = fields[1:]
                if rarefaction is None:
                    rarefaction = dict()
                    for sample in samples:
                        rarefaction[sample] = dict()
                        rarefaction[sample]['nb_asv'] = len([ i for i in biom.get_sample_obs(sample) if i >0 ])
                        rarefaction[sample]['nb_seq'] = biom.get_sample_count( sample )
                for sample in samples:
                    rarefaction[sample][rank] = list()
            else:
                if rarefaction_step_size is None:
                    rarefaction_step_size = int(fields[0])
                if rank not in rarefaction[sample]:
                    rarefaction[sample][rank] = list()
                for idx, sample in enumerate(samples):
                    if fields[idx+1] != "":
                        rarefaction[sample][rank].append( int(fields[idx+1]) )
        FH_rarefaction.close()
    del biom
    
    in_biom = BiomIO.from_json( input_biom )

    mode='Removed'
    if params.mask:
        mode = 'Masked'

    global_results = {
        'cluster': {
            'Kept': 0,
            'Modified': 0,
             mode: 0},
        'sequence' : {
            'Kept': 0,
            'Modified': 0,
             mode: 0,
        }
    }

    ### By sample and by filters
    samples_results = dict()
    filters_results = dict()

    ### track taxon lost
    # Number of taxon in each rank before filters
    taxon_lost = {'Blast' : {'all': dict() , 'multihit' : dict()}}
    for rank in range(len(params.taxonomic_ranks)):
        taxon_lost['Blast']["all"][rank] = list()
        taxon_lost['Blast']["multihit"][rank] = list()
    if in_biom.has_observation_metadata('rdp_taxonomy'):
        taxon_lost['RDP'] = {rank : list() for rank in range(len(params.taxonomic_ranks))}

    # store observation_name in dictionnary with all associated filter name (except Blast_taxonomy_changed)
    filters_intersections = dict()
    for filter in discards.keys():
        if filter != "Blast_taxonomy_changed" :
            for observation_name in discards[filter]:
                if not "Blast_taxonomy_changed" in discards or observation_name not in discards['Blast_taxonomy_changed']:
                    if observation_name not in filters_intersections:
                        filters_intersections[observation_name] = dict()
                    filters_intersections[observation_name][filter] = 1
    
    #initialized samples details results
    for sample in in_biom.get_samples_names():
        if not sample in samples_results:
            samples_results[sample] = {
                'initial': sum( 1 for x in in_biom.get_observations_by_sample(sample) ),
                'filtered': dict(),
                'kept': 0
            }
    # count number of observation filtered by a each combination of filters
    # store observation filtered by each criteria for each sample
    for observation_name in filters_intersections:
        # Removed intersection
        intersections_key = "--@@--".join(sorted( filters_intersections[observation_name].keys() ))
        if intersections_key not in filters_results:
            filters_results[intersections_key] = {
                'filters': list(filters_intersections[observation_name].keys()),
                'count': 0
            }
        filters_results[intersections_key]['count'] += 1

        # Filters by samples
        for sample in in_biom.get_samples_names():
            for filter in filters_intersections[observation_name]:
                if filter not in samples_results[sample]['filtered']:
                    samples_results[sample]['filtered'][filter] = 0
                if in_biom.get_count(observation_name,sample) > 0 :
                    samples_results[sample]['filtered'][filter] += 1

    # compute globale_results
    for observation_name in in_biom.get_observations_names():
        # ASV removed or taxonomy masked
        if observation_name in filters_intersections:
            global_results['cluster'][mode] += 1
            global_results['sequence'][mode] += in_biom.get_observation_count( observation_name )
        # blast taxonomy updated but not masked
        elif 'Blast_taxonomy_changed' in discards and observation_name in discards['Blast_taxonomy_changed']:
            global_results['cluster']['Modified'] += 1
            global_results['sequence']['Modified'] += in_biom.get_observation_count( observation_name )
            
            for sample in in_biom.get_samples_names():
                samples_results[sample]['kept'] += 1 if in_biom.get_count(observation_name,sample) > 0 else 0
        # keep unchanged
        else:
            global_results['cluster']['Kept'] += 1
            global_results['sequence']['Kept'] += in_biom.get_observation_count( observation_name )
            for sample in in_biom.get_samples_names():
                samples_results[sample]['kept'] += 1 if in_biom.get_count(observation_name,sample) > 0 else 0

        # track rdp taxon
        if 'RDP' in taxon_lost and in_biom.get_observation_metadata(observation_name)['rdp_taxonomy'] is not None:
            rdp_taxonomy = in_biom.get_observation_metadata(observation_name)['rdp_taxonomy']
            if issubclass(rdp_taxonomy.__class__,str):
                rdp_taxonomy = rdp_taxonomy.split(';')
            for i in range(len(params.taxonomic_ranks)):
                if ';'.join(rdp_taxonomy[:i+1]) not in taxon_lost['RDP'][i]:
                    taxon_lost['RDP'][i].append(';'.join(rdp_taxonomy[:i+1]))
        
        # track blast taxon
        if in_biom.get_observation_metadata(observation_name)['blast_affiliations'] is not None:
            for blast_affi in in_biom.get_observation_metadata(observation_name)['blast_affiliations'] :
                blast_taxonomy = blast_affi['taxonomy']
                if issubclass(blast_taxonomy.__class__,str):
                    blast_taxonomy = blast_taxonomy.split(';')
                for i in range(len(params.taxonomic_ranks)):
                    if ';'.join(blast_taxonomy[:i+1]) not in taxon_lost['Blast']['all'][i]:
                        taxon_lost['Blast']["all"][i].append(';'.join(blast_taxonomy[:i+1])) 
                    
            blast_taxonomy = in_biom.get_observation_metadata(observation_name)['blast_taxonomy']
            if issubclass(blast_taxonomy.__class__,str):
                blast_taxonomy = blast_taxonomy.split(';')
            
            for i in range(len(params.taxonomic_ranks)):
                if "Multi-affiliation" in ';'.join(blast_taxonomy[:i+1]) and ';'.join(blast_taxonomy[:i+1]) not in taxon_lost['Blast']["multihit"][i]:
                    taxon_lost['Blast']["multihit"][i].append(';'.join(blast_taxonomy[:i+1]))

    del in_biom

    # Remove taxon from taxon_lost if present in output biom
    out_biom = BiomIO.from_json( output_biom )
    for observation_name in out_biom.get_observations_names():
        # track lost rdp taxon
        if 'RDP' in taxon_lost and out_biom.get_observation_metadata(observation_name)['rdp_taxonomy'] is not None:
            rdp_taxonomy = out_biom.get_observation_metadata(observation_name)['rdp_taxonomy']
            if rdp_taxonomy:
                if issubclass(rdp_taxonomy.__class__,str):
                    rdp_taxonomy = rdp_taxonomy.split(';')
                for i in range(len(params.taxonomic_ranks)):
                    if ';'.join(rdp_taxonomy[:i+1]) in taxon_lost['RDP'][i]:
                        taxon_lost['RDP'][i].remove(';'.join(rdp_taxonomy[:i+1]))

        # track lost blast taxon
        if len(out_biom.get_observation_metadata(observation_name)['blast_affiliations'])>0:
            for blast_affi in out_biom.get_observation_metadata(observation_name)['blast_affiliations'] :
                blast_taxonomy = blast_affi['taxonomy']
                if issubclass(blast_taxonomy.__class__,str):
                    blast_taxonomy = blast_taxonomy.split(';')
                for i in range(len(params.taxonomic_ranks)):
                    if ';'.join(blast_taxonomy[:i+1]) in taxon_lost['Blast']['all'][i]:
                        taxon_lost['Blast']["all"][i].remove(';'.join(blast_taxonomy[:i+1])) 
                    
            blast_taxonomy = out_biom.get_observation_metadata(observation_name)['blast_taxonomy']
            if issubclass(blast_taxonomy.__class__,str):
                blast_taxonomy = blast_taxonomy.split(';')
            for i in range(len(params.taxonomic_ranks)):
                if "Multi-affiliation" in ';'.join(blast_taxonomy[:i+1]) and ';'.join(blast_taxonomy[:i+1]) in taxon_lost['Blast']["multihit"][i]:
                    taxon_lost['Blast']["multihit"][i].remove(';'.join(blast_taxonomy[:i+1]))

    del out_biom
    
    # store number instead of list of taxon lost
    for i in range(len(params.taxonomic_ranks)):
        taxon_lost['Blast'][params.taxonomic_ranks[i]] = [len(taxon_lost['Blast']['all'][i]) ,len(taxon_lost['Blast']['multihit'][i])]
        taxon_lost['Blast']['all'].pop(i)
        taxon_lost['Blast']['multihit'].pop(i)
        if 'RDP' in taxon_lost:
            taxon_lost['RDP'][params.taxonomic_ranks[i]] = len(taxon_lost['RDP'][i])
            taxon_lost['RDP'].pop(i)
    
    taxon_lost['Blast'].pop('all')
    taxon_lost['Blast'].pop('multihit')
    
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
        bootstrap_results = get_bootstrap_distrib( output_biom, args.bootstrap_tag, args.multiple_tag )
    # Get alignment metrics
    aln_results = get_alignment_distrib( output_biom, "perc_identity", "perc_query_coverage", "blast_affiliations" )

    # Write summary results
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "affiliation_filters_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###PROCESSED_FILTERS###" in line:
            line = line.replace( "###PROCESSED_FILTERS###", json.dumps([filter for filter in discards if filter != 'Blast_taxonomy_changed']) )
        if "###MODE###" in line:
            if "<h2 class=" in line:
                line = line.replace( "###MODE###", mode.lower() )    
            else:
                line = line.replace( "###MODE###", json.dumps(mode) )
        elif "###GLOBAL_RESULTS###" in line:
            line = line.replace( "###GLOBAL_RESULTS###", json.dumps(global_results) )
        elif "###GLOBAL_TAXON_LIST###" in line:
            line = line.replace( "###GLOBAL_TAXON_LIST###", json.dumps(params.taxonomic_ranks) )
        elif "###GLOBAL_TAXON_LOST###" in line:
            line = line.replace( "###GLOBAL_TAXON_LOST###", json.dumps(taxon_lost) )
        elif "###SAMPLES_RESULTS###" in line:
            line = line.replace( "###SAMPLES_RESULTS###", json.dumps(samples_results) )
        elif "###SAMPLES_NAMES###" in line:
            line = line.replace( "###SAMPLES_NAMES###", json.dumps(ordered_samples_names) )
        elif "###FILTERS_RESULTS###" in line:
            line = line.replace( "###FILTERS_RESULTS###", json.dumps(list(filters_results.values())) )
        elif "Draw a Venn to see which ASVs had been deleted by the filters chosen (Maximum 6 options): " in line and params.mask:
            line = "Draw a Venn to see which ASVs had its taxonomy masked by the filters chosen (Maximum 6 options): "
        elif "###TREE_DISTRIBUTION###" in line:
            line = line.replace( "###TREE_DISTRIBUTION###", json.dumps(newick_tree) )
        elif "###DATA_RAREFACTION###" in line:
            line = line.replace( "###DATA_RAREFACTION###", json.dumps(rarefaction) )
        elif "###RAREFACTION_STEP_SIZE###" in line:
            line = line.replace( "###RAREFACTION_STEP_SIZE###", json.dumps(rarefaction_step_size) )
        elif "###RAREFACTION_RANKS###" in line:
            line = line.replace( "###RAREFACTION_RANKS###", json.dumps(rarefaction_ranks) )
        elif "###ALIGNMENT_SCORES###" in line:
            line = line.replace( "###ALIGNMENT_SCORES###", json.dumps(aln_results) )
        elif "###BOOTSTRAP_SCORES###" in line:
            line = line.replace( "###BOOTSTRAP_SCORES###", json.dumps(bootstrap_results) )
        elif "###SAMPLES_DATA###" in line:
            line = line.replace( "###SAMPLES_DATA###", json.dumps(samples_results) )
        elif "###TAXONOMY_RANKS###" in line:
            line = line.replace( "###TAXONOMY_RANKS###", json.dumps(args.taxonomic_ranks) )
        elif "###FROGS_VERSION###" in line:
            line = line.replace( "###FROGS_VERSION###", "\""+str(__version__)+"\"" )
        elif "###FROGS_TOOL###" in line:
            line = line.replace( "###FROGS_TOOL###", "\""+ os.path.basename(__file__)+"\"" )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()



def process( args ):
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )

    impacted_biom = tmpFiles.add('impacted.biom')
    impacted_multihit = tmpFiles.add('impacted.multihit.tsv')
    impacted_biom2tsv_log= tmpFiles.add('impacted.biom2tsv.log')

    try:
        # parse biom, store impacted, write output biom by deleting ASV or masking taxonomies
        impacted_dict = filter_biom(args.input_biom,impacted_biom, args.output_biom, args)
        
        # write log
        Logger.static_write(args.log_file, "Identify ASV with :\n")
        for label in impacted_dict:
            if label != "Blast_taxonomy_changed" : 
                Logger.static_write(args.log_file, "\t- " + label + " : "+ str(len(impacted_dict[label])) + "\n")
        if 'Blast_taxonomy_changed' in impacted_dict:
            Logger.static_write(args.log_file,"\tadditionnaly, blast_taxonomy updated for " + str(len(impacted_dict['Blast_taxonomy_changed'])) + ' ASV(s)\n')

        # convert impacted biom in TSV
        if len(impacted_dict) > 0:
            BIOM_to_TSV(impacted_biom, args.input_fasta, args.impacted, args.impacted_multihit,impacted_biom2tsv_log).submit(args.log_file)
        else:
            BIOM_to_TSV(args.input_biom, args.input_fasta, args.impacted, args.impacted_multihit,impacted_biom2tsv_log, True).submit(args.log_file)

        # if params.delete, update fasta
        if args.delete:
            update_fasta_log = tmpFiles.add( "update_fasta_log.txt" )
            UpdateFasta( args.output_biom, args.input_fasta, args.output_fasta, update_fasta_log ).submit( args.log_file )
        
        # Affiliation stats

        tax_depth = [args.taxonomic_ranks.index(rank) for rank in args.taxonomic_ranks]
        rarefaction_ranks = args.taxonomic_ranks
        rarefaction_cmd = Rarefaction(args.output_biom, tmpFiles, "blast_taxonomy", tax_depth)
        rarefaction_cmd.submit( args.log_file )
        rarefaction_files = rarefaction_cmd.output_files
        # Put ASVs rarefaction file at the end , after species 
        rarefaction_files.append(rarefaction_files.pop(0))
        # Taxonomy tree
        tree_count_file = tmpFiles.add( "taxCount.enewick" )
        tree_ids_file = tmpFiles.add( "taxCount_ids.tsv" )
        TaxonomyTree(args.output_biom, "blast_taxonomy", tree_count_file, tree_ids_file).submit( args.log_file )
        
        # write summary
        write_summary( args.html, args.input_biom, args.output_biom, impacted_dict, tree_count_file, tree_ids_file, rarefaction_files, rarefaction_ranks, args )

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
    parser = argparse.ArgumentParser(description='Filters an abundance biom file on affiliations metrics')
    parser.add_argument('--version', action='version', version=__version__ )
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )
    
    parser.add_argument('--taxonomic-ranks', nargs='+', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels used in the metadata taxonomy. [Default: %(default)s]' )
    group_exclusion_taxonomy = parser.add_mutually_exclusive_group()
    group_exclusion_taxonomy.add_argument( '--taxonomy-tag', type=str, help='The metadata tag used in BIOM file to store the taxonomy. Use this parameter if the taxonomic affiliation has been processed by a software that adds only one affiliation or if you does not have a metadata with the consensus taxonomy (see "--tax-consensus-tag").Not allowed with --tax-consensus-tag. [Default: %(default)s]' )
    group_exclusion_taxonomy.add_argument( '--tax-consensus-tag', type=str, default="blast_taxonomy", help='The metadata tag used in BIOM file to store the consensus taxonomy. This parameter is used instead of "--taxonomy-tag" when you have several affiliations for each ASV. [Default: %(default)s]' )
    parser.add_argument( '--multiple-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the list of possible taxonomies. Use this parameter if the taxonomic affiliation has been processed by a software that adds several affiliation in the BIOM file (example: same score ambiguity). [Default: %(default)s]' )
    parser.add_argument( '--bootstrap-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the taxonomy bootstraps. [Default: %(default)s]' )
    parser.add_argument( '--identity-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the alignment identity. [Default: %(default)s]' )
    parser.add_argument( '--coverage-tag', type=str, default=None, help='The metadata tag used in BIOM file to store the alignment observation coverage. [Default: %(default)s]' )
    #     Filters behavior
    group_filter_bh = parser.add_argument_group( 'Filters behavior' )
    group_exclusion_filter_bh = group_filter_bh.add_mutually_exclusive_group()
    group_exclusion_filter_bh.add_argument('--mask', default=False, action='store_true', help="If affiliations do not respect one of the filter they are replaced by NA (mutually exclusive with --delete) [Default: %(default)s]")
    group_exclusion_filter_bh.add_argument('--delete', default=False, action='store_true', help="If affiliations do not respect one of the filter the entire ASV is deleted.(mutually exclusive with --mask) [Default: %(default)s]")
    #     Filters
    group_filter = parser.add_argument_group( 'Filters' )
    group_filter_blast_taxa = group_filter.add_mutually_exclusive_group()
    group_filter_blast_taxa.add_argument('--ignore-blast-taxa', type=str, nargs='*', help="Taxon list to masks/delete in Blast affiliations")
    group_filter_blast_taxa.add_argument('--keep-blast-taxa', type=str, nargs='*', help="Taxon list to keep in Blast affiliations. All others affiliations will be masks/delete.")
    group_filter.add_argument('--min-rdp-bootstrap', type=str, action=BootstrapParameter, metavar=("TAXONOMIC_LEVEL:MIN_BOOTSTRAP"), help="The minimal RDP bootstrap must be superior to this value (between 0 and 1)." )
    group_filter.add_argument('--min-blast-identity', type=ratioParameter, help="The number corresponding to the blast percentage identity (between 0 and 100)." )
    group_filter.add_argument('--min-blast-coverage', type=ratioParameter, help="The number corresponding to the blast percentage coverage (between 0 and 100)." )
    group_filter.add_argument('--max-blast-evalue', type=float, help="The number corresponding to the blast e value (between 0 and 1).")
    group_filter.add_argument('--min-blast-length', type=int, default=None, required=False, help="The number corresponding to the blast length." )
    #     Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-biom', required=True, help="The input biom file.")
    group_input.add_argument('--input-fasta', required=True, help="The input fasta file.")
    #     Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-biom', default="affiliation-filtered.biom", help="The Biom file output. [Default: %(default)s]")
    group_output.add_argument('--output-fasta', default="affiliation-filtered.fasta", help="The fasta output file. [Default: %(default)s]")
    group_output.add_argument('--html', default="summary.html", help="The HTML file containing the graphs. [Default: %(default)s]")
    group_output.add_argument('--impacted', default="impacted_clusters.tsv", help="The abundance file that summarizes all the clusters impacted (deleted or with affiliations masked). [Default: %(default)s]")
    group_output.add_argument('--impacted-multihit', default="impacted_clusters_multihit.tsv", help="The multihit TSV file associated with impacted ASV. [Default: %(default)s]")
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # keep quote around each taxon ignored when printing command line into logfile
    cmd=list()
    taxon = False
    for arg in sys.argv:
        if arg == "--ignore-blast-taxa" :
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
        raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : You must precise if you want to mask affiliations of delete ASV with --mask or --delete options.\n\n"))

    if args.min_rdp_bootstrap is None and args.min_blast_length is None and args.max_blast_evalue is None and args.min_blast_identity is None and args.min_blast_coverage is None:
        if args.ignore_blast_taxa is None and args.keep_blast_taxa is None:
            raise_exception(Exception("\n\n#ERROR : You need to specify at least on filtering criteria\n\n"))
        elif args.ignore_blast_taxa is not None and len(args.ignore_blast_taxa) == 0:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : At least one filter must be set to run " + os.path.basename(sys.argv[0]) + "\n\n"))
        elif args.keep_blast_taxa is not None and len(args.keep_blast_taxa) == 0:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : At least one filter must be set to run " + os.path.basename(sys.argv[0]) + "\n\n"))
    
    in_biom = BiomIO.from_json( args.input_biom )

    if args.min_rdp_bootstrap is not None:
        if not in_biom.has_observation_metadata("rdp_bootstrap"):
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'rdp_bootstrap'. You cannot use the parameter '--min-rdp-bootstrap' on this file.\n\n" ))
        elif not args.min_rdp_bootstrap["rank"] in args.taxonomic_ranks:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The taxonomic rank choosen in '--min-rdp-bootstrap' must be in '--taxonomic-ranks' (" + ";".join(args.taxonomic_ranks) + ").\n\n" ))
        else:
            nb_rank = 0
            for observation in in_biom.get_observations():
                if 'rdp_bootstrap' in observation['metadata'] and observation['metadata']['rdp_bootstrap'] is not None and len(observation['metadata']['rdp_bootstrap']) > 0:
                    nb_rank = len(observation['metadata']['rdp_bootstrap'])
                    break
            if nb_rank != len(args.taxonomic_ranks):
                raise_exception( argparse.ArgumentTypeError("\n\n#ERROR : RDP taxonomic metadata is defined on " + str(nb_rank) + " and you precise " + str(len(args.taxonomic_ranks)) + " ranks name (see --taxonomic-ranks)\n\n"))
    if not in_biom.has_observation_metadata("blast_affiliations"):
        if args.min_blast_length is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-length' on this file.\n\n" ))
        if args.max_blast_evalue is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-evalue' on this file.\n\n" ))
        if args.min_blast_identity is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--min-blast-identity' on this file.\n\n" ))
        if args.min_blast_coverage is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--max-blast-coverage' on this file.\n\n" ))
        if args.ignore_blast_taxa is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--ignore-blast-taxa' on this file.\n\n" ))
        if args.keep_blast_taxa is not None:
            raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : The BIOM input does not contain the metadata 'blast_affiliations'. You cannot use the parameter '--keep-blast-taxa' on this file.\n\n" ))     
    # Control non empty ignore-blast-taxa:
    if args.ignore_blast_taxa is not None:
        checkBlastTaxa(args.ignore_blast_taxa)
    if args.keep_blast_taxa is not None:
        checkBlastTaxa(args.keep_blast_taxa)       

    #for observation in in_biom.get_observations():
    #    taxonomy = observation['metadata']['blast_taxonomy']
    #    if taxonomy == None or len(taxonomy) == 0:
    #            print('\n\n#WARNING: you declare that taxonomies are defined on ' + str(len(args.taxonomic_ranks)) + ' ranks but your biom file contains taxonomy defined on ' + str(len(taxonomy)) + ', at least for ' + observation['id'] + '\n')
    #            print('Those clusters will be delete if --delete mode activated\n')
    #            break
    #del in_biom
    
    

    if args.delete and (not args.input_fasta or not args.output_fasta):
        raise_exception(Exception("\n\n#ERROR : In deletion mode, you must specify an input and output_fasta file\n\n"))

    # Process
    process( args )
