#!/usr/bin/env python3
#
# Copyright (C) 2022 INRAE
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

__author__ = ' Moussa Samb & Vincent Darbot & Geraldine Pascal - GENPHYSE '
__copyright__ = 'Copyright (C) 2022 INRAE'
__license__ = 'GNU General Public License'
__version__ = '4.0.1'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'dev'

import os
import sys
import json
import argparse
import pandas as pd

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

if os.getenv('GENE_HIERARCHY_FILE'):
   GENE_HIERARCHY_FILE=os.environ['GENE_HIERARCHY_FILE']
else:
   GENE_HIERARCHY_FILE=os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "frogsfunc_suppdata/gene_family_hierarchy.tsv"))

from frogsUtils import *
from frogsSequenceIO import *
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class MetagenomePipeline(Cmd):
	"""
	@summary: Per-sample metagenome functional profiles are generated based on the predicted functions for each study sequence.
	"""
	def __init__(self, in_biom, marker, function, max_nsti, min_reads, min_samples, strat_out, output_dir, log):
		"""
		@param in_biom: [str] Path to BIOM input file used in frogsfunc_placeseqs.
		@param marker: [str] Table of predicted marker gene copy numbers (frogsfunc_copynumbers output : frogsfunc_copynumbers_marker_nsti_predicted.tsv).
		@param function: [str] Table of predicted gene family copy numbers (frogsfunc_copynumbers output : frogsfunc_copynumbers_predicted_functions.tsv).
		@param max_nsti: [float] Sequences with NSTI values above this value will be excluded .
		@param min_reads: [int] Minimum number of reads across all samples for each input OTU.
		@param min_samples [int] Minimum number of samples that an OTU needs to be identfied within.
		@param strat_out: [boolean] if strat_out, output table stratified by sequences as well.
		@param function_abund: [str] Output file for metagenome predictions abundance.
		@param seqtab: [str] Output file with abundance normalized per marker copies number.
		@param weighted: [str] Output file with the mean of nsti value per sample.
		@param contrib: [str] Stratified output that reports contributions to community-wide abundances.
		"""
		opt = ' --strat_out ' if strat_out else ''

		Cmd.__init__(self,
				 'metagenome_pipeline.py ',
				 'Per-sample functional profiles prediction.',
				 " --input " +  in_biom + " --marker " + marker + " --function " + function + " --out_dir " + output_dir + " --max_nsti " + str(max_nsti) + " --min_reads " + str(min_reads) + " --min_samples " + str(min_samples) + opt + ' 2> ' + log,
				"--version")

	def get_version(self):
		 return "PICRUSt2 " + Cmd.get_version(self, 'stdout').split()[1].strip()

class ParseMetagenomePipeline(Cmd):
	"""
	@summary: Parse results of PICRUSt2 metageome_pipeline.py software to rerieve additional informations (i.g. databases functions links)
	"""
	def __init__(self, out_dir, out_abund, out_seqtab, out_weighted, strat_out, out_contrib, log):
		if strat_out:
			opt = " --input-contrib " + out_contrib 
		else:
			opt = ''
		Cmd.__init__( self,
					  'frogsFuncUtils.py',
					  'Parse metagenome_pipeline.py outputs.',
					  "parse-metagenome --input-dir " + out_dir + " --input-abund " + out_abund + " --input-seqtab " + out_seqtab + " --input-weighted " + out_weighted + opt + " 2>> " + log,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

class Biom2tsv(Cmd):
	"""
	@summary: Converts BIOM file to TSV file.
	"""
	def __init__(self, in_biom, out_tsv):
		Cmd.__init__( self,
					  'biom2tsv.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-file " + in_biom + " --output-file " + out_tsv + " --fields @observation_name @sample_count",
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

class Tsv2biom(Cmd):
	"""
	@summary: Create a temporary biom file that links every gene to samples abundances.
	This is necessary in order to display sunburst plots.
	"""
	def __init__(self, in_tsv, out_biom):

		Cmd.__init__( self,
					  'tsv_to_biom.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-tsv " + in_tsv + " --output-biom " + out_biom,
					  '--version' )

		self.in_tsv = in_tsv

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

	def parser(self, log_file):
		f_in = pd.read_csv(self.in_tsv, sep='\t')
		sum_col = f_in.pop("observation_sum")
		f_in.to_csv(self.in_tsv ,sep='\t' ,index=False)

class FormateAbundances(Cmd):
	"""
	@summary: Formate function abundances file in order to add function classifications and display sunbursts graphs.
	"""
	def __init__(self, in_abund, tmp_abund, hierarchy_file, log):

		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Formate function abundances file.',
			'formate-abundances --input-abundances ' + in_abund + ' --input-tmp-abundances ' + tmp_abund + ' --hierarchy-file ' + hierarchy_file + ' 2>> ' + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class TaxonomyTree(Cmd):
	"""
	@summary: Produces a tree with gene abundances by sample in extended newick format.
	"""
	def __init__(self, in_biom, taxonomy_tag, out_tree, out_ids):
		"""
		@param in_biom: [str] The processed BIOM path.
		@param taxonomy_tag: [str] The metadata title for the taxonomy in BIOM file.
		@param out_tree: [str] Path to the enewick output.
		@param out_ids: [str] Path to the IDs/samples output.
		"""
		Cmd.__init__( self,
					  'biomTools.py',
					  'Produces a taxonomy tree with counts by sample.',
					  'treeCount --input-file ' + in_biom + ' --taxonomy-key "' + taxonomy_tag + '" --output-enewick ' + out_tree + ' --output-samples ' + out_ids,
					  '--version' )

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def blast_threshold(in_biom, min_blast_identity, min_blast_coverage, out_biom):
	"""
	@summary: Removes sequences from biom file that does not pass the selected blast threshold.
	@param in biom: Biom file from frogsfunc_placeseqs step.
	@param blast_identity: Threshold of minimal blast % identity between the cluster sequence and its PICRUSt2 reference sequence.
	@param blast_coverage: Threshold of minimal blast % coverage between the cluster sequence and its PICRUSt2 reference sequence.
	@param out_biom: Output biom file without excluded sequences.
	"""
	biom = BiomIO.from_json(in_biom)
	discards = list()
	excluded_infos = dict()

	for observation in biom.get_observations():
		if min_blast_identity:
			ident = float(observation['metadata']['blast_picrust_ref_perc_identity'].split()[0]) / 100
			if ident < min_blast_identity:
				discards.append(observation['id'])
				excluded_infos[observation['id']] = dict()
				excluded_infos[observation['id']]['FROGS_taxonomy'] = str(';'.join(biom.get_observation_metadata(observation['id'])['blast_taxonomy']))
				excluded_infos[observation['id']]['PICRUSt2_taxonomy'] = biom.get_observation_metadata(observation['id'])['picrust2_affiliations']
				excluded_infos[observation['id']]['exclusion_paramater'] = "blast_identity"
				excluded_infos[observation['id']]['value_parameter'] = str(ident)
				continue
		if min_blast_coverage:
			cov = float(observation['metadata']['blast_picrust_ref_perc_query_coverage'].split()[0]) / 100
			if cov < min_blast_coverage:
				discards.append(observation['id'])
				excluded_infos[observation['id']] = dict()
				excluded_infos[observation['id']]['FROGS_taxonomy'] = str(';'.join(biom.get_observation_metadata(observation['id'])['blast_taxonomy']))
				excluded_infos[observation['id']]['PICRUSt2_taxonomy'] = biom.get_observation_metadata(observation['id'])['picrust2_affiliations']
				excluded_infos[observation['id']]['exclusion_paramater'] = "blast_coverage"
				excluded_infos[observation['id']]['value_parameter'] = str(cov)

	biom.remove_observations( discards )
	BiomIO.write( out_biom, biom )
	return excluded_infos

def excluded_sequence(in_biom, in_marker, out_seqtab, already_excluded, out_excluded_file):
	"""
	@summary: Returns the excluded sequence, that have a NSTI score above the NSTI threshold.
	@param in_biom: Biom file from frogsfunc_placeseqs step.
	@param in_marker: [str] Path to frogsfunc_copynumbers marker file to process.
	@param out_seqtab: [str] Path to frogsfunc_functions seqtab file to process.
	@output: The file of excluded sequence names.
	"""
	marker_file = open( in_marker ).readlines()[1:]
	seqtab_file = open( out_seqtab )
	biom = BiomIO.from_json(in_biom)
	excluded = open(out_excluded_file, "wt")
	cluster_to_nstis = dict()
	write_header = True
	no_excluded = True
	already_excluded_clusters = list()
	if len(already_excluded) > 0:
		excluded.write('\t'.join(['#Cluster','FROGS_taxonomy','PICRUSt2_taxonomy','exclusion_paramater','value_parameter'])+"\n")
		write_header = False
		no_excluded = False		
		for cluster, infos in already_excluded.items():
			already_excluded_clusters.append(cluster)
			excluded.write(cluster + "\t")
			for info in infos:
				excluded.write(already_excluded[cluster][info] + "\t")
			excluded.write('\n')

	clusters_in = list()
	for li in marker_file:
		cluster = li.strip().split('\t')[0]
		if cluster not in already_excluded_clusters:
			clusters_in.append(cluster)
			cluster_to_nstis[cluster] = li.strip().split('\t')[2]
	clusters_out = [ li.strip().split('\t')[0] for li in seqtab_file.readlines()[1:]]
	for cluster in clusters_in:
		if cluster not in clusters_out:
			no_excluded = False
			if write_header:
				excluded.write('\t'.join(['#Cluster','FROGS_taxonomy','PICRUSt2_taxonomy','exclusion_paramater','value_parameter'])+"\n")
				write_header = False
			excluded.write(cluster+"\t")
			try:
				excluded.write("\t".join([str(';'.join(biom.get_observation_metadata(cluster)['blast_taxonomy']))  ,str(biom.get_observation_metadata(cluster)['picrust2_affiliations']), "NSTI", cluster_to_nstis[cluster]])+"\n")
			except:
				excluded.write("\t".join(["unknown"  ,str(biom.get_observation_metadata(cluster)['picrust2_affiliations']), "NSTI", cluster_to_nstis[cluster]])+"\n")
	if no_excluded:
		excluded.write('#No excluded OTUs.\n')
	excluded.close()
	seqtab_file.close()

def write_summary(in_biom, function_file, nsti_file, excluded, tree_count_file, tree_ids_file, summary_file):
	"""
	@summary: Writes the process summary in one html file.
	@param in_biom: [str] path to the input BIOM file.
	@param function_file: [str] path to the gene abondancies fonction file.
	@param excluded: [str] The file of excluded sequence names.
	@param tree_count_file: [str] newick file of functions abundances per semple and per hierarchy.
	@param tree_ids_file: [str] file that link id to its sample.
	@param summary_file: [str] path to the output html file.
	"""
	# to summary OTUs number && abundances number
	summary_info = {
	   'nb_kept' : 0,
	   'nb_removed' : 0,
	   'abundance_kept' : 0,
	   'abundance_removed' : 0
	}
	number_otu_all = 0
	number_abundance_all = 0

	biom=BiomIO.from_json(in_biom)
	for otu in biom.get_observations_names():
		number_otu_all +=1
		number_abundance_all += biom.get_observation_count(otu)
	excluded_clusters = open( excluded ).readlines()
	if not excluded_clusters[0].startswith('#No excluded OTU'):
		#[1:] for skip header
		for otu in excluded_clusters[1:]:
			summary_info['nb_removed'] +=1
			summary_info['abundance_removed'] += biom.get_observation_count(otu.strip().split('\t')[0])

	summary_info['nb_kept'] = number_otu_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']

	samples_distrib = dict()
	FH_nsti = open(nsti_file).readlines()
	for li in FH_nsti:
		li = li.strip().split('\t')
		if li[0] in biom.get_samples_names():
			samples_distrib[li[0]] = {
			'mean_nsti' : round(float(li[1]), 3)
			}

	FH_tree_count = open( tree_count_file )
	newick_tree = FH_tree_count.readline()
	FH_tree_count.close()
	ordered_samples_names = list()
	FH_tree_ids = open( tree_ids_file )
	for line in FH_tree_ids:
		id, sample_name = line.strip().split( "\t", 1 )
		ordered_samples_names.append( sample_name )
	FH_tree_ids.close()

	# function abundances table
	infos_otus = list()
	details_categorys =["Function", "Description" ,"Observation_sum"]

	abund = open(function_file).readlines()
	#Header
	header = abund[0].strip().split('\t')
	for sample in header[4:]:
		details_categorys.append(sample)

	for li in abund[1:]:
		li = li.strip().split('\t')
		function = li[2]
		for i in range(len(li[3:])):
			sample = abund[0].strip().split('\t')[i+3]
			li[i+2] = round(float(li[i+3]),1)

		infos_otus.append({
			'name': li[2],
			'data': list(map(str,li[3:]))
			})

	# record details about removed OTU
	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "frogsfunc_functions_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
		elif "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_otus) )
		elif "###REMOVE_DATA###" in line:
			line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
		elif "###TAXONOMIC_RANKS###" in line:
			line = line.replace( "###TAXONOMIC_RANKS###", json.dumps(HIERARCHY_RANKS) )
		elif "###SAMPLES_NAMES###" in line:
			line = line.replace( "###SAMPLES_NAMES###", json.dumps(ordered_samples_names) )
		elif "###DATA_SAMPLE###" in line:
			line = line.replace( "###DATA_SAMPLE###", json.dumps(samples_distrib) )
		elif "###TREE_DISTRIBUTION###" in line:
			line = line.replace( "###TREE_DISTRIBUTION###", json.dumps(newick_tree) )
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
	parser = argparse.ArgumentParser( description='Per-sample functional profiles prediction.' )
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
	parser.add_argument('--strat-out', default=False, action='store_true', help='If activated, a new table is built. It will contain the abundances of each function of each OTU in each sample.')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input-biom', required=True, type=str, help='frogsfunc_placeseqs Biom output file (frogsfunc_placeseqs.biom).')
	group_input.add_argument('-f', '--input-function', required=True, type=str, help='Table of predicted gene family copy numbers (frogsfunc_copynumbers output, frogsfunc_copynumbers_predicted_functions.tsv).')
	group_input.add_argument('-m', '--input-marker', required=True, type=str, help='Table of predicted marker gene copy numbers (frogsfunc_copynumbers output, ex frogsfunc_copynumbers_marker.tsv).')
	group_input.add_argument('--max-nsti', type=float, default=2.0, help='Sequences with NSTI values above this value will be excluded (default: %(default)d).')
	group_input.add_argument('--min-blast-ident', type=float, default=None, help='Sequences with blast percentage identity against the PICRUSt2 closest ref above this value will be excluded (between 0 and 1)')
	group_input.add_argument('--min-blast-cov', type=float, default=None, help='Sequences with blast percentage coverage against the PICRUSt2 closest ref above this value will be excluded (between 0 and 1)')
	group_input.add_argument('--min-reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for each input OTU. OTUs below this cut-off will be counted as part of the \"RARE\" category in the stratified output. If you choose 1, none OTU will be grouped in “RARE” category.(default: %(default)d).')
	group_input.add_argument('--min-samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an OTU needs to be identfied within. OTUs below this cut-off will be counted as part of the \"RARE\" category in the stratified output.  If you choose 1, none OTU will be grouped in “RARE” category. (default: %(default)d).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-d', '--output-dir', default='frogsfunc_function_results', help='Output directory for functions predictions.')
	group_output.add_argument('--output-function-abund', default='frogsfunc_functions_unstrat.tsv', help='Output file for metagenome predictions abundance. (default: %(default)s).')
	group_output.add_argument('--output-seqtab', default='frogsfunc_functions_marker_norm.tsv', help='Output file with abundance normalized per marker copies number. (default: %(default)s).')
	group_output.add_argument('--output-weighted', default='frogsfunc_functions_weighted_nsti.tsv', help='Output file with the mean of nsti value per sample (format: TSV). [Default: %(default)s]' )
	group_output.add_argument('--output-contrib', default=None, help=' Stratified output that reports contributions to community-wide abundances (ex pred_metagenome_contrib.tsv).')
	group_output.add_argument('-e', '--excluded', default='frogsfunc_functions_excluded.txt', help='List of sequences with NSTI values above NSTI threshold ( --max_NSTI NSTI ).[Default: %(default)s]')
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--summary', default='frogsfunc_functions_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)
	tmp_files=TmpFiles(os.path.split(args.summary)[0])

	if args.strat_out and args.output_contrib is None:
		args.output_contrib = args.output_dir + '/frogsfunc_functions_strat.tsv'


	if not args.strat_out and args.output_contrib is not None:
		parser.error('--contrib FILENAME must be include with --strat_out flag')
	
	if args.min_blast_ident:
		if args.min_blast_ident < 0.0 or args.min_blast_ident > 1.0:
			parser.error('--min-blast-ident must be between 0.0 and 1.0.')
	if args.min_blast_cov:
		if args.min_blast_cov < 0.0 or args.min_blast_cov > 1.0:
			parser.error('--min-blast-cov must be between 0.0 and 1.0.')
	
	args.output_function_abund = args.output_dir + "/" + args.output_function_abund
	args.output_seqtab = args.output_dir + "/" + args.output_seqtab
	args.output_weighted = args.output_dir + "/" + args.output_weighted
	args.excluded = args.output_dir + "/" + args.excluded
	args.summary = args.output_dir + "/" + args.summary

	HIERARCHY_RANKS = ["Level1", "Level2", "Level3", "Gene"]
	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		excluded_infos = dict()
		if args.min_blast_ident or args.min_blast_cov:
			tmp_biom_blast_thresh = tmp_files.add( 'tmp_biom_blast_thresh' )
			tmp_excluded = tmp_files.add( 'tmp_excluded' )
			excluded_infos = blast_threshold(args.input_biom, args.min_blast_ident, args.min_blast_cov, tmp_biom_blast_thresh)

			tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
			Biom2tsv(tmp_biom_blast_thresh, tmp_biom_to_tsv).submit( args.log_file )
			biom_file = tmp_biom_blast_thresh
		else:
			#temp tsv file necessary for metagenome_pipeline.py
			tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
			Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )
			biom_file = args.input_biom

		tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )
		MetagenomePipeline(tmp_biom_to_tsv, args.input_marker, args.input_function, args.max_nsti, args.min_reads, args.min_samples, args.strat_out, args.output_dir, tmp_metag_pipeline).submit( args.log_file )
		
		tmp_parse = tmp_files.add( 'tmp_parse_metagenome.log' )
		ParseMetagenomePipeline(args.output_dir, args.output_function_abund, args.output_seqtab, args.output_weighted, args.strat_out, args.output_contrib, tmp_parse).submit( args.log_file)
		excluded_sequence(biom_file, args.input_marker, args.output_seqtab, excluded_infos, args.excluded)
		# Make a temporary functions abundances file to display sunbursts graphs.
		tmp_function_abund = tmp_files.add( "functions_unstrat.tmp")
		tmp_formate_abundances = tmp_files.add( 'tmp_formate_abundances.log' )
		FormateAbundances(args.output_function_abund, tmp_function_abund, GENE_HIERARCHY_FILE, tmp_formate_abundances).submit( args.log_file)
		

		tmp_biom = tmp_files.add( 'gene_abundances.biom' )
		Tsv2biom(tmp_function_abund, tmp_biom).submit( args.log_file)
		tree_count_file = tmp_files.add( "geneCount.enewick" )
		tree_ids_file = tmp_files.add( "geneCount_ids.tsv" )
		hierarchy_tag = "classification"
		TaxonomyTree(tmp_biom, hierarchy_tag, tree_count_file, tree_ids_file).submit( args.log_file )

		write_summary(args.input_biom, args.output_function_abund, args.output_weighted, args.excluded, tree_count_file, tree_ids_file, args.summary)
	finally:
		if not args.debug:
			tmp_files.deleteAll()
