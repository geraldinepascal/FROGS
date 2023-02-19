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

class HspFunction(Cmd):
	"""
	@summary: Predict number of marker copies (16S, 18S or ITS) for each cluster sequence (i.e OTU).
	"""
	def __init__(self, tree, marker_type, marker_file, function_table, functions, hsp_method, output_dir, output_file, nb_cpus, log_file):
		"""
		@param observed_marker_table: [str] Path to marker table file if marker studied is not 16S.
		@param in_tree: [str] Path to resulting tree file with inserted clusters sequences from frogsfunc_placeseqs.
		@param hsp_method: [str] HSP method to use.
		@param output: [str] PICRUSt2 marker output file.
		"""
		if marker_type != "16S":
			opt = ' --input-function-table ' + function_table
		elif function_table is None:
			opt = ' --functions ' + functions

		Cmd.__init__(self,
				 'launch_hsp.py',
				 'predict gene copy number per sequence.', 
				 ' function --input-tree ' + tree + ' --marker-type ' + marker_type + opt + ' --marker-file ' + marker_file + ' --hsp-method ' + hsp_method + ' --output-dir ' + output_dir + ' -o ' + output_file + ' --nb-cpus ' + str(nb_cpus) + '  --log-file ' + log_file,
				"--version")

		self.output = output_file
		self.log_file = log_file

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()


class MetagenomePipeline(Cmd):
	"""
	@summary: Per-sample metagenome functional profiles are generated based on the predicted functions for each study sequence.
	"""
	def __init__(self, in_biom, marker, function, max_nsti, min_reads, min_samples, strat_out, output_dir, log):
		"""
		@param in_biom: [str] Path to BIOM input file used in frogsfunc_placeseqs.
		@param marker: [str] Table of predicted marker gene copy numbers (frogsfunc_copynumbers output : frogsfunc_copynumbers_marker_nsti_predicted.tsv).
		@param function: [str] Table of predicted function copy numbers (frogsfunc_copynumbers output : frogsfunc_copynumbers_predicted_functions.tsv).
		@param max_nsti: [float] Sequences with NSTI values above this value will be excluded .
		@param min_reads: [int] Minimum number of reads across all samples for each input OTU.
		@param min_samples [int] Minimum number of samples that an OTU needs to be identfied within.
		@param strat_out: [boolean] if strat_out, output table stratified by sequences as well.
		@param function_abund: [str] Output file for function predictions abundance.
		@param otu_norm: [str] Output file with abundance normalized per marker copies number.
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
	def __init__(self, in_dir, out_abund, otu_norm_file , out_weighted, strat_out, out_contrib, log):
		opt = ''
		if strat_out:
			opt += " --output-contrib " + out_contrib
		Cmd.__init__( self,
					  'frogsFuncUtils.py',
					  'Parse metagenome_pipeline.py outputs.',
					  "parse-metagenome --input-dir " + in_dir + " --output-abund " + out_abund + " --output-seqtab " + otu_norm_file  + " --output-weighted " + out_weighted + opt + " 2>> " + log,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()


class RemoveSeqsBiomFasta(Cmd):
	'''
	@summary: Create new biom and fasta file without not insert sequences in tree.
	'''
	def __init__(self, in_fasta, in_biom, out_fasta, out_biom, excluded_file):
		'''
		@param in_fasta: [str] Path to fasta input file.
		@param in_biom: [str] Path to BIOM input file.
		@param out_fasta: [str] Path to fasta output file.
		@param out_biom: [str] Path to BIOM output file.
		@param excluded_file: [str] Path to not insert sequences file (Cluster ID in the first column).
		'''
		Cmd.__init__(self,
			'remove_seqs_biom_fasta.py',
			'remove not insert sequences in tree from fasta and biom file.',
			'--input-biom ' + in_biom + ' --input-fasta ' + in_fasta + ' --excluded-sequences ' + excluded_file + ' --output-biom ' + out_biom + " --output-fasta " + out_fasta,
			'--version')

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
	def __init__(self, in_abund, tmp_unstrat, hierarchy_file, log):

		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Formate function abundances file.',
			'formate-abundances --input-abundances ' + in_abund +  ' --input-tmp-unstrat ' + tmp_unstrat + ' --hierarchy-file ' + hierarchy_file + ' 2>> ' + log,
			'--version')

	def get_version(self):
		return Cmd.get_version(self, 'stdout').strip()

class GenerateSunburst(Cmd):
	"""
	@summary: Generate sunburst input files for html graphics
	"""
	def __init__(self, in_abund, tmp_sunburst, log):

		Cmd.__init__(self,
			'frogsFuncUtils.py',
			'Generate sunburst input files.',
			'generate-sunburst --input-abundances ' + in_abund + ' --input-tmp-sunburst ' + tmp_sunburst + ' 2>> ' + log,
			'--version')

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

def otus_filter(in_biom, nsti_file, min_blast_identity, min_blast_coverage, max_nsti, excluded_file):
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
	nstis = dict()
	FH_excluded = open(excluded_file, 'wt')

	with open(nsti_file) as FH_nsti:
		next(FH_nsti)
		for li in FH_nsti:
			nstis[li.split('\t')[0]] = float(li.strip().split('\t')[-1])
	
	for observation in biom.get_observations():
		cov = float(observation['metadata']['blast_picrust_ref_perc_query_coverage'].split()[0]) / 100
		ident = float(observation['metadata']['blast_picrust_ref_perc_identity'].split()[0]) / 100
		nsti = nstis[observation['id']]
		exclusion_paramater = str()
		value_paramater = str()
		if min_blast_identity and ident < min_blast_identity:
			exclusion_paramater = "min_blast_identity"
			value_paramater = "identity = " + str(ident)
		if min_blast_coverage and cov < min_blast_coverage:
			if exclusion_paramater == "":
				exclusion_paramater = "min_blast_coverage"
				value_paramater = "coverage = " + str(cov)
			else:
				exclusion_paramater += ",min_blast_coverage"
				value_paramater += ",coverage = " + str(cov)
		if max_nsti and nsti > max_nsti:
			if exclusion_paramater == "":
				exclusion_paramater = "max_nsti"
				value_paramater = "nsti = " + str(nsti)
			else:
				exclusion_paramater += ",max_nsti"
				value_paramater += ",nsti = " + str(nsti)	
		if exclusion_paramater != "":
			discards.append(observation['id'])
			excluded_infos[observation['id']] = dict()
			excluded_infos[observation['id']]['FROGS_taxonomy'] = str(';'.join(biom.get_observation_metadata(observation['id'])['blast_taxonomy']))
			excluded_infos[observation['id']]['PICRUSt2_taxonomy'] = biom.get_observation_metadata(observation['id'])['picrust2_affiliations']
			excluded_infos[observation['id']]['exclusion_paramater'] = exclusion_paramater
			excluded_infos[observation['id']]['value_parameter'] = value_paramater
	if len(excluded_infos) > 0:
		FH_excluded.write('\t'.join(['#Cluster','FROGS_taxonomy','PICRUSt2_taxonomy','exclusion_paramater','value_parameter'])+"\n")
		for excluded in excluded_infos:
			FH_excluded.write('\t'.join([excluded, excluded_infos[excluded]['FROGS_taxonomy'], excluded_infos[excluded]['PICRUSt2_taxonomy'], excluded_infos[excluded]['exclusion_paramater'], excluded_infos[excluded]['value_parameter']]) + "\n")
	else:
		FH_excluded.write('#No excluded OTUs.\n')
	FH_excluded.close()
	return excluded_infos

def check_nsti_threshold(max_nsti, in_biom):
	'''
	Test if the NSTI threshold specified by the user is less than the minimum NSTI on the dataset.
	'''
	biom = BiomIO.from_json(in_biom)
	min_nsti = None
	for observation in biom.get_observations():
		if biom.get_observation_metadata(observation['id'])['NSTI']:
			cur_nsti = float(biom.get_observation_metadata(observation['id'])['NSTI'])
			if min_nsti == None:
				min_nsti = cur_nsti
			elif cur_nsti < min_nsti:
				min_nsti = cur_nsti

	if args.max_nsti < min_nsti:
		return raise_exception( Exception( "\n\n#ERROR : --max-nsti " + str(max_nsti) + " threshold will remove all clusters.\n\n" ))

def check_basename_files(arg_name, file_path):
	'''
	Test if output file name specified by the user only contains the file name, without directory.
	'''
	if not os.path.basename(file_path) == file_path:
		return raise_exception( Exception( "\n\n#ERROR : --" + arg_name.replace('_','-') + \
		" should only contain a filename, without directory (You specified " + arg_value +" ). Please use --output-dir to specify the output directory.\n\n"))

def count_nb_obs_per_ranks(in_biom):
	'''
	rank_to_obs associates each taxonomic level rank to its observations.
	@Output: List with number of different observations per taxonomic rank.
	'''
	rank_to_obs = {
		"0" : [],
		"1" : [],
		"2" : [],
		"3" : [],
		"4" : [],
		"5" : [],
		"6" : [],
	}
	no_info_obs = ["Multi-affiliation", "unknown species"]
	biom = BiomIO.from_json(in_biom)

	for observation in biom.get_observations():
		if biom.has_metadata("blast_taxonomy"):
			taxo_hiera = biom.get_observation_metadata(observation['id'])["blast_taxonomy"]
			for i in range(len(taxo_hiera)):
				if taxo_hiera[i] not in rank_to_obs[str(i)] and taxo_hiera[i] not in no_info_obs:
					rank_to_obs[str(i)].append(taxo_hiera[i])
	return [ len(rank_to_obs[str(i)]) for i in rank_to_obs ]

def write_summary(in_biom, function_file, nsti_file, excluded, tree_count_file, tree_ids_file, out_biom, summary_file):
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
	# Construct star_plot about number of different taxonomics ranks retrieved before and after different thresholds.
	in_taxo_ranks, out_taxo_ranks = count_nb_obs_per_ranks(args.input_biom), count_nb_obs_per_ranks(args.output_biom)
	starplot_series = {
		'before_series' : in_taxo_ranks,
		'after_series' : out_taxo_ranks
	}

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
		elif "###STARPLOT_SERIES###" in line:
			line = line.replace( "###STARPLOT_SERIES###", json.dumps(starplot_series) )
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
	parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
	parser.add_argument('--strat-out', default=False, action='store_true', help='If activated, a new table is built. It will contain the abundances of each function of each OTU in each sample.')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input-biom', required=True, type=str, help='frogsfunc_placeseqs Biom output file (frogsfunc_placeseqs.biom).')
	group_input.add_argument('-i', '--input-fasta', required=True, help='frogsfunc_placeseqs Fasta output file (frogsfunc_placeseqs.fasta).')
	group_input.add_argument('-t', '--input-tree', required=True, type=str, help='frogsfunc_placeseqs output tree in newick format containing both studied sequences (i.e. ASVs or OTUs) and reference sequences.')
	group_input.add_argument('-m', '--input-marker', required=True, type=str, help='Table of predicted marker gene copy numbers (frogsfunc_placeseqs output : frogsfunc_marker.tsv).')
	group_input.add_argument('--marker-type', required=True, choices=['16S','ITS','18S'], help='Marker gene to be analyzed.')
	
	group_input_16S = parser.add_argument_group( '16S ' )
	group_input_16S.add_argument('-f', '--functions', default=["EC"], nargs='+', choices=['EC', 'KO', 'COG', 'PFAM', 'TIGRFAM','PHENO'], help="Specifies which function databases should be used (%(default)s). EC is used by default because necessary for frogsfunc_pathways. At least EC or KO is required. To run the command with several functions, separate the functions with spaces (ex: -i EC PFAM).")
	group_input_other = parser.add_argument_group( 'ITS and 18S ' )
	group_input_other.add_argument('--input-function-table',help="The path to input functions table describing directly observed functions, in tab-delimited format.(ex $PICRUSt2_PATH/default_files/fungi/ec_ITS_counts.txt.gz). Required.")
    
	group_input.add_argument('--hsp-method', default='mp', choices=['mp', 'emp_prob', 'pic', 'scp', 'subtree_average'], help='HSP method to use. mp: predict discrete traits using max parsimony. emp_prob: predict discrete traits based on empirical state probabilities across tips. subtree_average: predict continuous traits using subtree averaging. pic: predict continuous traits with phylogentic independent contrast. scp: reconstruct continuous traits using squared-change parsimony (default: %(default)s).')
	group_input.add_argument('--max-nsti', type=float, default=2.0, help='Sequences with NSTI values above this value will be excluded (default: %(default)d).')
	group_input.add_argument('--min-blast-ident', type=float, default=None, help='Sequences with blast percentage identity against the PICRUSt2 closest ref above this value will be excluded (between 0 and 1)')
	group_input.add_argument('--min-blast-cov', type=float, default=None, help='Sequences with blast percentage coverage against the PICRUSt2 closest ref above this value will be excluded (between 0 and 1)')
	group_input.add_argument('--min-reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for each input OTU. OTUs below this cut-off will be counted as part of the \"RARE\" category in the stratified output. If you choose 1, none OTU will be grouped in “RARE” category.(default: %(default)d).')
	group_input.add_argument('--min-samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an OTU needs to be identfied within. OTUs below this cut-off will be counted as part of the \"RARE\" category in the stratified output.  If you choose 1, none OTU will be grouped in “RARE” category. (default: %(default)d).')
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('-d', '--output-dir', default='frogsfunc_function_results', help='Output directory for function predictions.')
	group_output.add_argument('--output-function', default="frogsfunc_copynumbers_functions.tsv", type=str, help='Output table with predicted function abundances per studied sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.[Default: %(default)s]')
	group_output.add_argument('--output-function-abund', default='frogsfunc_functions_unstrat.tsv', help='Output file for function prediction abundances. (default: %(default)s).')
	group_output.add_argument('--output-otu-norm', default='frogsfunc_functions_marker_norm.tsv', help='Output file with otu abundances normalized by marker copies number. (default: %(default)s).')
	group_output.add_argument('--output-weighted', default='frogsfunc_functions_weighted_nsti.tsv', help='Output file with the mean of nsti value per sample (format: TSV). [Default: %(default)s]' )
	group_output.add_argument('--output-contrib', default=None, help=' Stratified output that reports otu contributions to community-wide function abundances (ex pred_function_otu_contrib.tsv).')
	group_output.add_argument('--output-biom', default='frogsfunc_function.biom', help='Biom file without excluded OTUs (NSTI, blast perc identity or blast perc coverage thresholds). (format: BIOM) [Default: %(default)s]')
	group_output.add_argument('--output-fasta', default='frogsfunc_function.fasta', help='Fasta file without excluded OTUs (NSTI, blast perc identity or blast perc coverage thresholds). (format: FASTA). [Default: %(default)s]')
	group_output.add_argument('-e', '--output-excluded', default='frogsfunc_functions_excluded.txt', help='List of OTUs with NSTI values above NSTI threshold ( --max_NSTI NSTI ).[Default: %(default)s]')
	group_output.add_argument('-l', '--log-file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-s', '--summary', default='frogsfunc_functions_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )
	args = parser.parse_args()
	prevent_shell_injections(args)
	args_dict = vars(args)
	for arg_name, arg_value in args_dict.items():
		if arg_name.startswith('output') and arg_name != "output_dir" and arg_value is not None:
			check_basename_files(arg_name, arg_value)
	if not os.path.exists(args.output_dir):
		os.mkdir(args.output_dir)
			
	### Check inputs
	# Check for 16S input
	if args.marker_type == "16S" and (not 'EC' in args.functions and not 'KO' in args.functions):
			parser.error("\n\n#ERROR : --input-functions : 'EC' and/or 'KO' must be at least indicated (others functions are optionnal)")
	# Check for ITS or 18S input
	if args.marker_type in ["ITS", "18S"]:
		if args.input_function_table is None:
			parser.error("\n\n#ERROR : --input-function-table required when studied marker is not 16S!\n\n")

	if not args.strat_out and args.output_contrib is not None:
		parser.error('--strat_out flag must be include with --output-contrib')
	
	if args.min_blast_ident:
		if args.min_blast_ident < 0.0 or args.min_blast_ident > 1.0:
			parser.error('--min-blast-ident must be between 0.0 and 1.0.')
	if args.min_blast_cov:
		if args.min_blast_cov < 0.0 or args.min_blast_cov > 1.0:
			parser.error('--min-blast-cov must be between 0.0 and 1.0.')
	###
	### Output paths

	args.output_otu_norm = args.output_dir + "/" + args.output_otu_norm
	args.output_weighted = args.output_dir + "/" + args.output_weighted
	args.output_fasta = args.output_dir + "/" + args.output_fasta
	args.output_biom = args.output_dir + "/" + args.output_biom
	args.output_excluded = args.output_dir + "/" + args.output_excluded
	args.summary = args.output_dir + "/" + args.summary
	if args.strat_out:
		if args.output_contrib is None:
			args.output_contrib = "frogsfunc_functions_strat.tsv"
		args.output_contrib = args.output_dir + "/" + args.output_contrib

	tmp_files=TmpFiles(os.path.split(args.output_otu_norm)[0])
	tmp_files_picrust =  TmpFiles(os.path.split(args.output_otu_norm)[0])

	HIERARCHY_RANKS = ["Level1", "Level2", "Level3", "Function_id"]
	try:
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

		check_nsti_threshold(args.max_nsti, args.input_biom)
		excluded_infos = dict()
		if args.min_blast_ident or args.min_blast_cov or args.max_nsti:
			tmp_biom_blast_thresh = tmp_files.add( 'tmp_biom_blast_thresh' )
			tmp_excluded = tmp_files.add( 'tmp_excluded' )
## ecrire ligne loger.static write sur l'exclusion des parametres en question
			excluded_infos = otus_filter(args.input_biom, args.input_marker, args.min_blast_ident, args.min_blast_cov, args.max_nsti, args.output_excluded)

			RemoveSeqsBiomFasta(args.input_fasta, args.input_biom, args.output_fasta, args.output_biom, args.output_excluded).submit(args.log_file)
			tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
			Biom2tsv(args.output_biom, tmp_biom_to_tsv).submit( args.log_file )

		else:
			#temp tsv file necessary for metagenome_pipeline.py
			tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
			Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )

		functions = " ".join(args.functions)
		tmp_hsp_function = tmp_files.add( 'tmp_hsp_function.log' )
		args.output_function = args.output_dir + "/" + args.output_function
		HspFunction(args.input_tree, args.marker_type, args.input_marker, args.input_function_table, functions, args.hsp_method, args.output_dir, args.output_function, args.nb_cpus, tmp_hsp_function).submit(args.log_file)
		FH_in = open(tmp_hsp_function)
		for line in FH_in:
			if line.startswith('## Software :'):
				tool_version =  line.strip().replace("##Software :", "##Software : PICRUSt2 ")
		FH_in.close()
		Logger.static_write(args.log_file, tool_version + "\n\n")

		function_outputs = [function + "_copynumbers_predicted.tsv" for function in args.functions]
		for function_file in function_outputs:
			database = function_file.split('_')[0]
			tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )
			function_file = args.output_dir + "/" + function_file
			##
			tmp_files_picrust =  TmpFiles(os.path.dirname(function_file), prefix="")
			tmp_seqtab = tmp_files_picrust.add('seqtab_norm.tsv.gz')
			tmp_weighted = tmp_files_picrust.add('weighted_nsti.tsv.gz')
			tmp_unstrat = tmp_files_picrust.add('pred_metagenome_unstrat.tsv.gz')
			output_strat_abund = None
			if args.strat_out:
				strat_basename_ext = os.path.basename(args.output_contrib)
				strat_basename = os.path.splitext(strat_basename_ext)[0]
				ext = os.path.splitext(strat_basename_ext)[1]
				output_strat_abund = args.output_dir + "/" + strat_basename + "_" + database + ext
				tmp_strat = tmp_files_picrust.add('pred_metagenome_contrib.tsv.gz')
			##
			MetagenomePipeline(tmp_biom_to_tsv, args.input_marker, function_file, args.max_nsti, args.min_reads, args.min_samples, args.strat_out, args.output_dir, tmp_metag_pipeline).submit( args.log_file )
			function_basename_ext = os.path.basename(args.output_function_abund)
			function_basename = os.path.splitext(function_basename_ext)[0]
			ext = os.path.splitext(function_basename_ext)[1]
			output_function_abund = args.output_dir + "/" + function_basename + "_" + database + ext
			tmp_parse = tmp_files.add( 'tmp_parse_metagenome.log' )
			ParseMetagenomePipeline(args.output_dir, output_function_abund, args.output_otu_norm, args.output_weighted, args.strat_out, output_strat_abund, tmp_parse).submit( args.log_file)
				
			tmp_function_unstrat = tmp_files.add( "functions_unstrat.tmp")
			tmp_formate_abundances = tmp_files.add( 'tmp_formate_abundances.log' )
			FormateAbundances(output_function_abund, tmp_function_unstrat, GENE_HIERARCHY_FILE, tmp_formate_abundances).submit( args.log_file)		
			
			# Launch one time for EC or KO.
			to_run = True
			if (database == "EC" or database == "KO") and to_run:
				# Make a temporary functions abundances file to display sunbursts graphs.
				tmp_function_sunburst = tmp_files.add( "functions_unstrat_sunburst.tmp")
				function_file_sunburst = output_function_abund
				tmp_sunburst_log = tmp_files.add( 'tmp_generate_sunburst.log' )
				GenerateSunburst(output_function_abund, tmp_function_sunburst, tmp_sunburst_log).submit( args.log_file)
				function_file_sunburst = tmp_function_sunburst
				to_run = False

		tmp_biom = tmp_files.add( 'gene_abundances.biom' )
		Tsv2biom(tmp_function_sunburst, tmp_biom).submit( args.log_file)
		tree_count_file = tmp_files.add( "geneCount.enewick" )
		tree_ids_file = tmp_files.add( "geneCount_ids.tsv" )
		hierarchy_tag = "classification"
		TaxonomyTree(tmp_biom, hierarchy_tag, tree_count_file, tree_ids_file).submit( args.log_file )

		write_summary(args.input_biom, function_file_sunburst, args.output_weighted, args.output_excluded, tree_count_file, tree_ids_file, args.output_biom, args.summary)
	finally:
		if not args.debug:
			tmp_files.deleteAll()
			tmp_files_picrust.deleteAll()
