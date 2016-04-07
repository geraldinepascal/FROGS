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
__version__ = '0.10.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'beta'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import Biom, BiomIO
from frogsNode import Node


####################################################################################################################
#
# Tree count
#
####################################################################################################################
def task_treeCount( args ):
    tree, ordered_samples = get_tree_with_count( args.input_file, (args.output_samples is not None), args.taxonomy_key )

    # Writes the extended newick
    FH_enewick = open( args.output_enewick, "w" )
    FH_enewick.write( tree.to_extended_newick() )
    FH_enewick.close()

    # Writes the link between samples names and samples IDs
    if args.output_samples is not None:
        FH_samples_ids = open( args.output_samples, "w" )
        for id, sample_name in enumerate( ordered_samples ):
            FH_samples_ids.write( str(id) + "\t" + str(sample_name) + "\n" )
        FH_samples_ids.close()

def get_tree_with_count( input_biom, compress=False, taxonomy_key="taxonomy" ):
    """
    @summary: Returns the tree of taxa and their counts by sample from BIOM.
    @param input_biom: [str] Path to the BIOM file processed.
    @param compress: [bool] if true the samples names are replaced by samples index.
    @param taxonomy_key: [str] The metadata title for the taxonomy in biom.
    @return: [list] The tree generated and the ordered list of samples names (usefull to retrieve name by index if you use compress).
    """
    ordered_samples_names = list()
    tree = Node("root")
    biom = BiomIO.from_json( input_biom )
    for sample_name in biom.get_samples_names():
        ordered_samples_names.append( sample_name )
        sample_id = None if not compress else (len(ordered_samples_names)-1)
        update_tree_for_sample( biom, tree, sample_name, taxonomy_key, sample_id )
    return tree, ordered_samples_names

def update_tree_for_sample( biom, tree, sample_name, taxonomy_key, sample_id=None ):
    """
    @summary: Updates a tree with sample taxa (nodes and counts).
    @param biom: [Biom] The Biom object.
    @param tree: [Node] The root node of the tree to update.
    @param sample_name: [str] The sample name to process.
    @param taxonomy_key: [str] The metadata title for the taxonomy in biom.
    @param sample_id: [str] The sample id to replace the sample name in tree.
    """
    sample_key = sample_name if sample_id is None else str(sample_id)
    for observation in biom.get_observations_by_sample( sample_name ):
        current_node = tree
        if observation["metadata"].has_key(taxonomy_key) and observation["metadata"][taxonomy_key] is not None:
            # Get taxonomy
            taxonomy = biom.get_observation_taxonomy( observation["id"], taxonomy_key )
            # Add taxon in tree
            for taxon in taxonomy:
                if not current_node.has_child( taxon ):
                    current_node.add_child( Node(taxon) )
                current_node = current_node.get_child( taxon )
            # Add sample count in node
            if not current_node.metadata.has_key(sample_key):
                current_node.metadata[sample_key] = 0
            current_node.metadata[sample_key] += biom.get_count( observation["id"], sample_name )
    return tree


####################################################################################################################
#
# Sampling
#
####################################################################################################################
def task_sampling( args ):
    if args.nb_sampled is None and args.sampled_ratio is None:
        raise Exception('--nb-sampled or --sampled-ratio must be provided.')
    sampling_by_sample( args.input_file, args.output_file, args.nb_sampled, args.sampled_ratio )


def sampling_by_sample( input_biom, output_biom, nb_sampled=None, sampled_ratio=None ):
    """
    @summary: Writes a BIOM after a random sampling in each sample.
    @param input_biom: [str] Path to the processed BIOM.
    @param output_biom: [str] Path to outputed BIOM.
    @param nb_sampled: [int] Number of sampled sequences by sample.
    @param sampled_ratio: [float] Ratio of sampled sequences by sample.
    @note: nb_sampled and sampled_ratio are mutually exclusive.
    """
    initial_biom = BiomIO.from_json( input_biom )
    new_biom = Biom(
                    matrix_type="sparse",
                    generated_by="Sampling " + (str(nb_sampled) if nb_sampled is not None else str(sampled_ratio) + "%" ) + " elements by sample from " + input_biom
    )
    observations_already_added = dict()
    for sample_name in initial_biom.get_samples_names():
        new_biom.add_sample( sample_name, initial_biom.get_sample_metadata(sample_name) )
        sample_seq = initial_biom.get_sample_count(sample_name)
        sample_nb_sampled = nb_sampled
        if nb_sampled is None:
            sample_nb_sampled = int(sample_seq * sampled_ratio)
        if sample_seq < nb_sampled:
            raise Exception( str(sample_nb_sampled) + " sequences cannot be sampled in sample '" + str(sample_name) + "'. It only contains " + str(sample_seq) + " sequences." )
        else:
            for current_nb_iter in range(sample_nb_sampled):
                # Take an observation in initial BIOM
                selected_observation = initial_biom.random_obs_by_sample(sample_name)
                selected_observation_id = selected_observation['id']
                initial_biom.subtract_count( selected_observation_id, sample_name, 1 )
                # Put in new BIOM
                if not observations_already_added.has_key(selected_observation_id):
                    new_biom.add_observation( selected_observation_id, initial_biom.get_observation_metadata(selected_observation_id) )
                    observations_already_added[selected_observation_id] = True
                new_biom.add_count( selected_observation_id, sample_name, 1 )
    BiomIO.write( output_biom, new_biom )


####################################################################################################################
#
# Rarefaction
#
####################################################################################################################
def task_rarefaction( args ):
    rarefaction_data = rarefaction( args.input_file, args.step_size, args.ranks, args.taxonomy_key )
    for current_rank in args.ranks:
        output_file = args.output_file_pattern.replace('##RANK##', str(current_rank))
        write_output( output_file, rarefaction_data[current_rank], args.step_size )


def rarefaction( input_biom, interval=10000, ranks=None, taxonomy_key="taxonomy" ):
    """
    @summary: Returns the rarefaction by ranks by samples.
    @param input_biom: [str] Path to the biom file processed.
    @param interval: [int] Size of first sampling.
    @param ranks: [list] The rank(s) level for the diversity.
                   Example :
                     Sampled set :
                       Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Sphingomonas
                       Bacteria; Proteobacteria; Gammaproteobacteria; Vibrionales; Vibrionaceae; Vibrio; Vibrio halioticoli
                       Bacteria; Proteobacteria; Gammaproteobacteria; Legionellales; Coxiellaceae; Coxiella; Ornithodoros moubata symbiont A
                       Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Burkholderiaceae; Limnobacter; Limnobacter thiooxidans
                     Result for this set
                       With rank 1 or 2 : 1 group
                       With rank 3 : 3 different groups
                       With rank 4 or 5 or 6 : 4 different groups
    @param taxonomy_key : [str] The metadata title for the taxonomy in the input.
    @return: [dict] By ranks by samples the list of differents taxa for each steps.
              Example :
                 {
                   1: {
                        "sampleA" : [10, 20, 22, 23, 24, 25, 25, 25 ],
                        "sampleB" : [15, 25, 28, 30, 32, 34, 35, 36, 37, 37, 37, 37]
                      }
                 }
    @warning: The taxa with name starting with unknown used as complete new name 'unknown'.
    """
    sample_rarefaction = dict()
    biom = BiomIO.from_json( input_biom )
    for current_rank in ranks:
        sample_rarefaction[current_rank] = dict()
    for sample in biom.get_samples_names():
        taxa = dict()
        for current_rank in ranks:
            sample_rarefaction[current_rank][sample] = list()
            taxa[current_rank] = dict()
        sample_count = biom.get_sample_count( sample )
        expected_nb_iter = sample_count/interval
        for current_nb_iter in range(expected_nb_iter):
            selected_observations = biom.random_obs_extract_by_sample(sample, interval)
            for current_selected in selected_observations:
                taxonomy = list()
                if current_selected['observation']["metadata"].has_key(taxonomy_key) and current_selected['observation']["metadata"][taxonomy_key] is not None:
                    taxonomy = biom.get_observation_taxonomy( current_selected['observation']["id"], taxonomy_key )
                for idx, taxon in enumerate(taxonomy):
                    if taxon.lower().startswith("unknown"):
                        taxonomy[idx] = "unknown"
                while len(taxonomy) < max(ranks):
                    taxonomy.append("unknown")
                for current_rank in ranks:
                    taxonomy_str = (';'.join(taxonomy[0:current_rank+1])).lower()
                    taxa[current_rank][taxonomy_str] = True
            for current_rank in ranks:
                sample_rarefaction[current_rank][sample].append( str(len(taxa[current_rank])) )
    return sample_rarefaction


def write_output( output_path, rarefaction_data, interval, joiner="\t" ):
    """
    @summary : Write rarefaction data in file.
    @param output_path : [str] Path to the output file.
    @param rarefaction_data : [dict] By sample the list of differents taxa for each steps (see function rarefaction).
    @param interval : [int] Size of first sampling.
    @param joiner : [str] The output fields separator.
    @note : Example of one output file
              #Nb_sampled    SampleA    SampleB
              100    8    4
              200    10    5
              300    12    5
              400    12    5
              500    12
    """
    FH_out = open( output_path, "w" )
    FH_out.write( "#Nb_sampled" + joiner + joiner.join(sorted(rarefaction_data.keys())) + "\n" )

    stop = False
    current_nb_iter = 0
    while not stop:
        stop = True
        line = str( (current_nb_iter + 1) * interval )
        nb_selected = (current_nb_iter + 1) * interval
        for sample in sorted(rarefaction_data.keys()):
            if len(rarefaction_data[sample]) >= (current_nb_iter + 1):
                stop = False
                line += joiner + str(rarefaction_data[sample][current_nb_iter])
            else:
                line += joiner
        if not stop:
            FH_out.write( line + "\n" )
        current_nb_iter += 1


####################################################################################################################
#
# Observations depth
#
####################################################################################################################
def task_observationDepth( args ):
    observations_depth( args.input_file, args.output_file )


def observations_depth( input_biom, output_depth ):
    """
    @summary : Write the depths of the observation in file.
    @param input_biom : [str] path to the biom file processed.
    @param output_depth : [str] path to the output file.
    @note : Example of one output file
                #Depth<TAB>Nb_Observ_concerned<TAB>Prct_Observ_concerned
                1<TAB>65<TAB>65.000
                2<TAB>30<TAB>30.000
                3<TAB>0<TAB>0.000
                4<TAB>5<TAB>5.000
    """
    obs_depth = list()
    nb_observ = 0
    # Process depth calculation
    biom = BiomIO.from_json( input_biom )
    for observation_id, observation_count in biom.get_observations_counts():
        while len(obs_depth) <= observation_count:
            obs_depth.append(0)
        obs_depth[observation_count] += 1
        if observation_count != 0:
            nb_observ += 1
    del biom
    # Write output
    out_fh = open( output_depth, 'w' )
    out_fh.write( "#Depth\tNb_Observ_concerned\tPrct_Observ_concerned\n" )
    for depth in range(1, len(obs_depth)):
        prct = (float(obs_depth[depth])/ nb_observ)*100
        out_fh.write( str(depth) + "\t" + str(obs_depth[depth]) + "\t" + ("%.3f" % prct) + "\n" )
    out_fh.close()


####################################################################################################################
#
# Hierarchical classification
#
####################################################################################################################
def task_hclassification( args ):
    samples_hclassification( args.input_file, args.output_file, args.distance_method, args.linkage_method )


def to_newick( node, id_2_name ):
    """
    @return: [str] The newick representation of the tree rooted by the node.
    @param node: [scipy.cluster.hierarchy.ClusterNode] The node use as root for newick tree.
    @param id_2_name: [dict] Names of sample by node ID.
    """
    return "(" + _to_newick_part(node, id_2_name) + ");"

def _to_newick_part( node, id_2_name ):
    """
    @return: [str] The newick representation of the node and her descendant.
    @param node: [scipy.cluster.hierarchy.ClusterNode] The node to convert in newick.
    @param id_2_name: [dict] Names of sample by node ID.
    """
    childrens = list()
    if node.left or node.right:
        if node.left: childrens.append( _to_newick_part(node.left, id_2_name) )
        if node.right: childrens.append( _to_newick_part(node.right, id_2_name) )
        return  "(" + ",".join( childrens ) + "):" + ("%.3f" % float(node.dist))
    else:
        return id_2_name[node.id]

def samples_hclassification( input_biom, output_newick, distance_method, linkage_method ):
    """
    @summary : Process and write an hierarchical classification from Biom.
    @param input_biom : [str] Path to the BIOM file to process.
    @param output_newick : [str] Path to the newick output file.
    @param distance_method : [str] Used distance method for classify.
    @param linkage_method : [str] Used linkage method for classify.
    """
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import linkage, dendrogram
    import scipy.cluster.hierarchy
    data_array = list()
    samples_names = list()

    # Normalisation on count by sample
    biom = BiomIO.from_json( input_biom )
    for col_idx, current_sample in enumerate(biom.columns):
        samples_names.append( current_sample['id'] )
        sum_on_sample = biom.data.get_col_sum( col_idx )
        OTUs_norm = list()
        for row_idx in range(len(biom.rows)):
            OTUs_norm.append( biom.data.nb_at(row_idx, col_idx)/float(sum_on_sample) )
        data_array.append( OTUs_norm )
    del biom

    if len(samples_names) == 1 :
        # Write newick
        out_fh = open( output_newick, "w" )
        out_fh.write( "(" + samples_names[0] + ");\n" )
        out_fh.close()
    else:
        # Computing the distance and linkage
        data_dist = pdist( data_array, distance_method )
        data_link = linkage( data_dist, linkage_method )
        # Write newick
        scipy_hc_tree = scipy.cluster.hierarchy.to_tree( data_link , rd=False )
        id_2_name = dict( zip(range(len(samples_names)), samples_names) )
        out_fh = open( output_newick, "w" )
        out_fh.write( to_newick(scipy_hc_tree, id_2_name) + "\n" )
        out_fh.close()


####################################################################################################################
#
# BIOM to TSV
#
####################################################################################################################
def task_biom2tsv( args ):
    biom_to_tsv( args.input_file, args.output_file, args.fields, args.list_separator )


def biom_to_tsv( input_biom, output_tsv, fields, list_separator ):
    """
    @summary: Convert BIOM file to TSV file.
    @param input_biom: [str] Path to the BIOM file.
    @param output_tsv: [str] Path to the output file (format : TSV).
    @param fields: [list] Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count'. The others columns must be metadata title.
    @param list_separator: [str] Separator for complex metadata.
    """
    biom = BiomIO.from_json( input_biom )
    out_fh = open( output_tsv, "w" )
    # Header
    line = list()
    for current_field in fields:
        if current_field == '@observation_name':
            line.append( "observation_name" )
        elif current_field == '@sample_count':
            line.append( "\t".join(biom.get_samples_names()) )
        elif current_field == '@observation_sum':
            line.append( "observation_sum" )
        else: #metadata
            line.append( str(current_field) )
    out_fh.write( "#" + "\t".join(line) + "\n" )
    # Data
    for idx, count_by_sample in enumerate(biom.to_count()):
        observation = biom.rows[idx]
        line = list()
        for current_field in fields:
            if current_field == '@observation_name':
                line.append( str(observation['id']) )
            elif current_field == '@sample_count':
                line.append( "\t".join(map(str, count_by_sample)) )
            elif current_field == '@observation_sum':
                line.append( str(sum(count_by_sample)) )
            else: #metadata
                if issubclass(observation['metadata'][current_field].__class__, list):
                    line.append( list_separator.join(observation['metadata'][current_field]) )
                else:
                    line.append( str(observation['metadata'][current_field]) )
        out_fh.write( "\t".join(line) + "\n" )
    out_fh.close()


####################################################################################################################
#
# Argparse types
#
####################################################################################################################
def strict_positive_int(value):
    """
    @summary: Raises an exception if the argument value is < 1.
    """
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("The minimum value is 1.")
    return value


####################################################################################################################
#
# Main
#
####################################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description='biomTools is a program package designed for working with Biom files.' )
    parser.add_argument( "--version", action='version', version=__version__ )
    subparsers = parser.add_subparsers()

    # Sampling parameters
    parser_sampling = subparsers.add_parser('sampling', help='Random sampling in each sample.', usage='''biomTools.py sampling [-h] -i INPUT_FILE -o OUTPUT_FILE
                             -n NB_SAMPLED | -r SAMPLED_RATIO''')
    parser_sampling.add_argument( '-i', '--input-file', required=True, type=str, help='BIOM file processed.' )
    parser_sampling.add_argument( '-o', '--output-file', required=True, type=str, help='Sampling results in BIOM format.' )
    group_exclusion_sampling = parser_sampling.add_mutually_exclusive_group()
    group_exclusion_sampling.add_argument( '-n', '--nb-sampled', type=int, help='Number of sampled sequences by sample.' )
    group_exclusion_sampling.add_argument( '-r', '--sampled-ratio', type=float, help='Ratio of sampled sequences by sample.' )
    parser_sampling.set_defaults(func=task_sampling)

    # Rarefaction parameters
    parser_rarefaction = subparsers.add_parser('rarefaction', help='Process data for rarefaction curve by sample.')
    parser_rarefaction.add_argument( '-i', '--input-file', required=True, type=str, help='BIOM file processed.' )
    parser_rarefaction.add_argument( '-o', '--output-file-pattern', required=True, type=str, help='Rarefaction file(s) pattern with tag "##RANK##". Example: "/tmp/rarefaction_##RANK##.tsv".' )
    parser_rarefaction.add_argument( '-s', '--step-size', type=strict_positive_int, default=10000, help='Additional number of sampled sequences by round of sampling.' )
    parser_rarefaction.add_argument( '-r', '--ranks', nargs='+', required=True, type=int, default=None, help='The taxonomy depth used to evaluate diversity.' )
    parser_rarefaction.add_argument( '-k', '--taxonomy-key', type=str, default="taxonomy", help='The metadata title for the taxonomy in your BIOM file. Example : "rdp_taxonomy"' )
    parser_rarefaction.set_defaults(func=task_rarefaction)

    # Hierarchical classification parameters
    parser_hclassification = subparsers.add_parser('hclassification', help='Process data for hierarchical classification dendrogram.')
    parser_hclassification.add_argument( '-i', '--input-file', required=True, type=str, help='BIOM file processed.' )
    parser_hclassification.add_argument( '-o', '--output-file', required=True, type=str, help='Hierarchical classification in NEWICK format.' )
    parser_hclassification.add_argument( '-d', '--distance-method', type=str, default="euclidean", help='Used distance method for classify (example : euclidean).' )
    parser_hclassification.add_argument( '-l', '--linkage-method', type=str, default="average", help='used linkage method for classify (example : centroid).' )
    parser_hclassification.set_defaults(func=task_hclassification)

    # Observation depth parameters
    parser_observationDepth = subparsers.add_parser('obsdepth', help='Process observation depth.')
    parser_observationDepth.add_argument( '-i', '--input-file', required=True, type=str, help='BIOM file processed.' )
    parser_observationDepth.add_argument( '-o', '--output-file', required=True, type=str, help='Depth file (Depth<TAB>Nb_Observ<TAB>Prct_Observ).' )
    parser_observationDepth.set_defaults(func=task_observationDepth)

    # Biom 2 tsv parameters
    parser_biom2tsv = subparsers.add_parser('biom2tsv', help='Convert BIOM file to TSV file.')
    parser_biom2tsv.add_argument( '-i', '--input-file', required=True, help='BIOM file processed.' )
    parser_biom2tsv.add_argument( '-o', '--output-file', required=True, help='Path to the output file (format : TSV).')
    parser_biom2tsv.add_argument( '-f', '--fields', default=['@observation_name', '@observation_sum', '@sample_count'], nargs='+', help="Columns and their order in output. Special columns : '@observation_name', '@observation_sum', '@sample_count'. The others columns must be metadata titles.")
    parser_biom2tsv.add_argument( '-s', '--list-separator', default=';', help='Separator for complex metadata.')
    parser_biom2tsv.set_defaults(func=task_biom2tsv)

    # Tree count parameters
    parser_treeCount = subparsers.add_parser('treeCount', help='Produces a taxonomy tree with counts by sample in extended newick format.')
    parser_treeCount.add_argument( '-i', '--input-file', required=True, help='BIOM file processed.' )
    parser_treeCount.add_argument( '-e', '--output-enewick', required=True, help='Path to the output file (format : enewick).')
    parser_treeCount.add_argument( '-s', '--output-samples', type=str, help="Path to the output file with link between samples names and ids (format : TSV). If this option is used the samples names in enewick are replaced by ids to reduce the file weight.")
    parser_treeCount.add_argument( '-k', '--taxonomy-key', type=str, default="taxonomy", help='The metadata title for the taxonomy in your BIOM file. Example : "rdp_taxonomy"' )
    parser_treeCount.set_defaults(func=task_treeCount)

    # Parse parameters and call process
    args = parser.parse_args()
    args.func(args)
