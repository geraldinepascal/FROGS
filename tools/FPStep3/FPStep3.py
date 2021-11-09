#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard & Vincent Darbot & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import os
import re
import sys
import gzip
import glob
import json
import shutil
import argparse
import pandas as pd
from tempfile import gettempdir

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
GENE_HIERARCHY_FILE = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "default_files/gene_family_hierarchy.tsv"))

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
	def __init__(self, in_biom, marker, function, max_nsti, min_reads, min_samples, strat_out, function_abund, seqtab, weighted, contrib, log):
		"""
		@param: in_biom: Path to BIOM input file used in FPStep1
		"""
		opt = ' --strat_out ' if strat_out else ''

		Cmd.__init__(self,
				 'metagenome_pipeline.py ',
				 'Per-sample functional profiles prediction.', 
				 " --input " +  in_biom + " --marker " + marker + " --function " + function + " --out_dir ./ --max_nsti " + str(max_nsti) + " --min_reads " + str(min_reads) + " --min_samples " + str(min_samples) + opt + ' 2> ' + log,
				"--version")

		self.abund = function_abund
		self.seqtab = seqtab
		self.weighted = weighted
		self.strat = strat_out
		self.contrib = contrib

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').split()[1].strip() 

	def parser(self, log_file):
		START_GENBANK_LINK = "https://www.genome.jp/dbget-bin/www_bget?"
		START_COG_LINK = "https://www.ncbi.nlm.nih.gov/research/cog/cog/"
		START_PFAM_LINK = "https://pfam.xfam.org/family/"
		START_TIGR_LINK = "https://0-www-ncbi-nlm-nih-gov.linyanti.ub.bw/genome/annotation_prok/evidence/"
		f_in = gzip.open('pred_metagenome_unstrat.tsv.gz', 'rt').readlines()
		f_out = open(self.abund, 'wt')
		header = f_in[0].strip().split('\t')
		header.insert(0,'db_link')
		f_out.write("\t".join(header)+"\n")
		for li in f_in[1:]:
			li = li.strip().split('\t')
			function = li[0]
			if "COG" in function:
				li.insert(0,START_COG_LINK + function )
			if "PF" in function:
				li.insert(0,START_PFAM_LINK + function )
			if "TIGR" in function:
				li.insert(0,START_TIGR_LINK + function )
			elif re.search('K[0-9]{5}',function) or "EC:" in function:
				li.insert(0,START_GENBANK_LINK + function )
			f_out.write("\t".join(li)+"\n")
		os.remove('pred_metagenome_unstrat.tsv.gz')
		with gzip.open('seqtab_norm.tsv.gz', 'rb') as f_in:
			with open(self.seqtab, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('seqtab_norm.tsv.gz')
		with gzip.open('weighted_nsti.tsv.gz', 'rb') as f_in:
			with open(self.weighted, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			os.remove('weighted_nsti.tsv.gz')
		if self.strat:
			with gzip.open('pred_metagenome_contrib.tsv.gz', 'rb') as f_in:
				with open(self.contrib, 'wb') as f_out:
					shutil.copyfileobj(f_in, f_out)
				os.remove('pred_metagenome_contrib.tsv.gz')

class Biom2tsv(Cmd):
	"""
	@summary: Converts BIOM file to TSV file.
	"""
	def __init__(self, in_biom, out_tsv):

		Cmd.__init__( self,
					  'biom2tsv.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-file " + in_biom + " --output-file " + out_tsv,
					  '--version' )

	def get_version(self):
		 return Cmd.get_version(self, 'stdout').strip()

class Tsv2biom(Cmd):
	"""
	@summary: In order to creates a temporary biom file that links every gene to samples abundances.
	This is necessary in order to display sunburst plots.

	"""
	def __init__(self, in_tsv, out_biom):

		Cmd.__init__( self,
					  'tsv_to_biom.py',
					  'Converts a BIOM file in TSV file.',
					  "--input-tsv " + in_tsv + " --output-biom " + out_biom,
					  '--version' )

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
		# Cmd
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

def excluded_sequence(in_biom, in_marker, out_seqtab, excluded):
	"""
	@summary: Returns the excluded sequence, that have a NSTI score above the NSTI threshold.
	@param in_biom: Biom file.
	@param in_marker: [str] Path to FPStep2 marker file to process.
	@param out_seqtab: [str] Path to FPStep3 seqtab file to process.
	@output: The file of excluded sequence names.
	"""
	marker_file = open( in_marker )
	seqtab_file = open( out_seqtab )
	biom = BiomIO.from_json(in_biom)
	excluded = open(excluded, "wt")
	clusters_in = [ li.strip().split('\t')[0] for li in marker_file.readlines()[1:]]
	clusters_out = [ li.strip().split('\t')[0] for li in seqtab_file.readlines()[1:]]
	no_excluded = True
	write_header = True
	for cluster in clusters_in:
		if cluster not in clusters_out:
			no_excluded = False
			if write_header:
				excluded.write('\t'.join(['Cluster','FROGS_taxonomy','Picrust2_taxonomy'])+"\n")
				write_header = False
			excluded.write(cluster+"\t")
			excluded.write("\t".join([str(';'.join(biom.get_observation_metadata(cluster)['blast_taxonomy']))  ,str(biom.get_observation_metadata(cluster)['picrust2_affiliations'])])+"\n")
	if no_excluded:
		excluded.write('#No excluded OTUs.\n')
	excluded.close()
	marker_file.close()
	seqtab_file.close()

def formate_abundances_file(strat_file, gene_hierarchy_file, hierarchy_tag = "hierarchy"):
	"""
	@summary: Formate FPSTep3 output in order to create a biom file of pathways abundances.
	@param strat_file: FPStep3 output of gene abundances prediction (FPStep3_pred_metagenome_unstrat.tsv)
	@param gene_hierarchy_file: reference file that links every gene ID to its hierarchy levels.
	@param tmp_tsv: temporary tsv output of abundances per samples.
	"""
	id_to_hierarchy = {}
	path_fi = open(gene_hierarchy_file).readlines()
	for li in path_fi:
		li = li.strip().split('\t')
		id_to_hierarchy[li[-1]] = ";".join(li)

	df = pd.read_csv(strat_file,sep='\t')
	df.rename(columns = {'function':'observation_name'}, inplace = True)
	headers = ['observation_name', 'db_link']
	for column in df:
		if column not in headers:
			df[column] = df[column].round(0).astype(int)

	df.to_csv(strat_file ,sep='\t' ,index=False)
	tmp = open(strat_file + ".tmp", 'wt')
	FH_in = open(strat_file).readlines()
	header = FH_in[0].strip().split('\t')
	header.insert(0, hierarchy_tag)
	tmp.write("\t".join(header)+"\n")
	for li in FH_in[1:]:
		li = li.strip().split('\t')
		if li[1] in id_to_hierarchy:
			li.insert(0,id_to_hierarchy[li[1]])
			tmp.write("\t".join(li)+"\n")
	tmp.close()
	os.rename(strat_file +'.tmp', strat_file)
	return hierarchy_tag

def write_summary(in_biom, strat_file, excluded, tree_count_file, tree_ids_file, summary_file):
	"""
	@summary: Writes the process summary in one html file.
	@param in_biom: [str] path to the input BIOM file.
	@param strat_file: [str] path to the gene abondancies fonction file.
	@param excluded: [str] The file of excluded sequence names.
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
	if not excluded_clusters[0].startswith('#'):
		#[1:] for skip header
		for otu in excluded_clusters[1:]:
			summary_info['nb_removed'] +=1
			summary_info['abundance_removed'] += biom.get_observation_count(otu.strip().split('\t')[0])

	summary_info['nb_kept'] = number_otu_all - summary_info['nb_removed']
	summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']

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

	abund = open(strat_file)
	for li in abund:
		if "observation_name" in li:
			li = li.strip().split('\t')
			for sample in li[4:]:
				details_categorys.append(sample)
			break

	for li in abund:
		li = li.strip().split('\t')
		function = li[2]
		for i in range(len(li[3:])):
			li[i+2] = round(float(li[i+3]),1)

		infos_otus.append({
			'name': li[2],
			'data': list(map(str,li[3:]))
			})
	abund.close()
	# record details about removed OTU

	FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep3_tpl.html") )
	FH_summary_out = open( summary_file, "wt" )

	for line in FH_summary_tpl:
		if "###DETECTION_CATEGORIES###" in line:
			line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(details_categorys) )
		elif "###DETECTION_DATA###" in line:
			line = line.replace( "###DETECTION_DATA###", json.dumps(infos_otus) )
		elif "###REMOVE_DATA###" in line:
			line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
		elif "###TAXONOMIC_RANKS###" in line:
			line = line.replace( "###TAXONOMIC_RANKS###", json.dumps(args.hierarchy_ranks) )
		elif "###SAMPLES_NAMES###" in line:
			line = line.replace( "###SAMPLES_NAMES###", json.dumps(ordered_samples_names) )
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
	parser.add_argument('--strat_out', default=False, action='store_true', help='Output table stratified by sequences as well. By default this will be in \"contributional\" format ''(i.e. long-format) unless the \"--wide_table\" ''option is set. The startified outfile is named ''\"metagenome_contrib.tsv.gz\" when in long-format.')
	# Inputs
	group_input = parser.add_argument_group( 'Inputs' )
	group_input.add_argument('-b', '--input_biom', required=True, type=str, help='FPStep1 sequence abundances output file (FPStep1.biom).')
	group_input.add_argument('-f', '--function', required=True, type=str, help='Table of predicted gene family copy numbers ''(FPStep2 output, ex FPStep2_all_predicted.tsv).')
	group_input.add_argument('-m', '--marker', required=True, type=str, help='Table of predicted marker gene copy numbers ''(FPStep2 output, ex FPStep2_marker_nsti_predicted.tsv.')
	group_input.add_argument('--max_nsti', type=float, default=2.0, help='Sequences with NSTI values above this value will ' 'be excluded (default: %(default)d).')
	group_input.add_argument('--min_reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for ''each input ASV. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	group_input.add_argument('--min_samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an ASV needs to be ''identfied within. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')
	group_input.add_argument('--hierarchy_ranks', nargs='*', default=["Level1", "Level2", "Level3", "Gene"], help='The ordered ranks levels used in the metadata hierarchy pathways. [Default: %(default)s]' )
	#Outputs
	group_output = parser.add_argument_group( 'Outputs')
	group_output.add_argument('--function_abund', default='FPStep3_pred_metagenome_unstrat.tsv', help='Output file for metagenome predictions abundance. (default: %(default)s).')
	group_output.add_argument('--seqtab', default='FPStep3_seqtab_norm.tsv', help='This output file will contain abundance normalized. (default: %(default)s).')
	group_output.add_argument('--weighted', default='FPStep3_weighted_nsti.tsv', help='This output file will contain the nsti value per sample (format: TSV). [Default: %(default)s]' )
	group_output.add_argument('--contrib', default=None, help=' Stratified output that represents contributions to community-wide abundances (ex pred_metagenome_contrib.tsv)')
	group_output.add_argument('-e', '--excluded', default='FPStep3_excluded.txt', help='List of sequences with NSTI values above NSTI threshold ( --max_NSTI NSTI ).')
	group_output.add_argument('-l', '--log_file', default=sys.stdout, help='List of commands executed.')
	group_output.add_argument('-t', '--html', default='FPStep3_summary.html', help="Path to store resulting html file. [Default: %(default)s]" )	
	args = parser.parse_args()
	prevent_shell_injections(args)

	if args.strat_out and args.contrib is None:
		args.contrib = 'FPStep3_pred_metagenome_contrib.tsv'

	if not args.strat_out and args.contrib is not None:
		parser.error('--contrib FILENAME must be include with --strat_out flag')

	tmp_files=TmpFiles(os.path.split(args.marker)[0])
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		#temp tsv file necessary for metagenome_pipeline.py
		tmp_biom_to_tsv = tmp_files.add( 'tmp_biom_to_tsv' )
		Biom2tsv(args.input_biom, tmp_biom_to_tsv).submit( args.log_file )

		tmp_metag_pipeline = tmp_files.add( 'tmp_metagenome_pipeline.log' )	
		MetagenomePipeline(tmp_biom_to_tsv, args.marker, args.function, args.max_nsti, args.min_reads, args.min_samples, args.strat_out, args.function_abund, args.seqtab, args.weighted, args.contrib, tmp_metag_pipeline).submit( args.log_file )
		
		excluded_sequence(args.input_biom, args.marker, args.seqtab, args.excluded)

		hierarchy_tag = formate_abundances_file(args.function_abund, GENE_HIERARCHY_FILE)
		tmp_biom = tmp_files.add( 'gene_abundances.biom' )
		Tsv2biom(args.function_abund, tmp_biom).submit( args.log_file)
		tree_count_file = tmp_files.add( "geneCount.enewick" )
		tree_ids_file = tmp_files.add( "geneCount_ids.tsv" )
		TaxonomyTree(tmp_biom, hierarchy_tag, tree_count_file, tree_ids_file).submit( args.log_file )

		write_summary(args.input_biom, args.function_abund, args.excluded, tree_count_file, tree_ids_file, args.html)
	finally:
		if not args.debug:
			tmp_files.deleteAll()