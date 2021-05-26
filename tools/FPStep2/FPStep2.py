#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

#Import
import os
import sys
import argparse
import json
import re
from numpy import median
from collections import OrderedDict
import pandas as pd
#import argparse
#from picrust2.wrap_hsp import castor_hsp_workflow
#from picrust2.util import make_output_dir_for_file, check_files_exist
#from picrust2.default import default_tables

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPAT
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
# Default table PATH
ITS = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/ITS_counts.txt.gz")
EC = os.path.join(os.path.dirname(os.__file__), "site-packages/picrust2/default_files/fungi/ec_ITS_counts.txt.gz")


#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class hsp(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: hsp.py program  predict number copie of gene family for each OTU.
	@summary: use 16S, EC and/or KO 
	@see: https://github.com/picrust/picrust2/wiki
	"""
	
	def __init__(self, in_trait, tree, category, output, log):
	#def __init__(self, hsp16Met, categorie, in_trait, tree, output, n, stdout):	
	   #os.system("hsp.py "+ hsp16Met  +" -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output) +" --observed_trait_table "+ str(observed_trait_table) + " -n " )
		#-n: Calculate NSTI and add to output file.
        if category == "16S":
            input_marker = " -i 16S "

        elif category == "ITS"


		Cmd.__init__(self,
				 'hsp.py',
				 'predict gene copy 16S', 
				  " -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output) +" -n" +' 2> ' + log,
				"--version") 

		
	def get_version(self):
		"""
		@summary: Returns the program version number.
		@return: [str] Version number if this is possible, otherwise this method return 'unknown'.
		"""

		return Cmd.get_version(self, 'stdout').split()[1].strip()
		

#####################################################   
class hspITS(Cmd):
	"""
	@summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
	@summary: hsp.py program  predict number copie of gene family for each OTU.
	@summary: use EC 
	@see: https://github.com/picrust/picrust2/wiki
	@commande: hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 1
	"""
	def __init__(self, hspITSMet,categorie, tree, output, observed_trait_table, n, stdout):
		#os.system("hsp.py "+ hsp_EC_Met  +" -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output))
		#os.system("hsp.py "+ hspITSMet+ " -t "+ str(tree) +" -o "+ str(output) +" --observed_trait_table "+ str(observed) + " -n ")
		ITS = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "hsp/picrust2/default_files/fungi/ITS_counts.txt.gz"))

	   
		Cmd.__init__(self,
				 'hsp.py',
				 'predict gene copy fungi', 
				  "" + hspITSMet  +" -t "+ str(tree) +" -o "+ str(output) + " --observed_trait_table "+ str(observed_trait_table) +" -n" +' 2> ' + stdout,
				"--version") 
	
	def get_version(self): 
		"""
		@summary: Returns the program version number.
		@return: [str] Version number if this is possible, otherwise this method return 'unknown'.
		"""

		return Cmd.get_version(self, 'stdout').split()[1].strip() 



##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_summary( summary_file, results_chimera ):
	"""
	@summary: Writes the summary of results.
	@param summary_file: [str] The output file.
	@param results_chimera: [str] Path to the input chimera step summary.
	"""
	# Get data
	detection_categories = ["Sequence", "Otu", "NSTI"]
	detection_data = []
	log_fh = open(results_chimera, "r")
	line = log_fh.readline()

	for line in log_fh:
		detection_data.append(line.strip().split("\t"))


	#detection_categories = ["Kept nb", "Kept abundance", "Removed nb", "Removed abundance", "Abundance of the most abundant removed", "Detected nb", "Detected abundance", "Abundance of the most abundant detected"]
	# detection_data = list()
	# remove_data = dict()

	# # Parse results chimera
	# in_remove_metrics = True
	# in_detection_metrics = False
	# section_first_line = True
	# log_fh = open(results_chimera)
	# for line in log_fh:
	#	 line = line.strip()
	#	 if line.startswith('##Metrics by sample'):
	#		 remove_metrics = False
	#		 in_detection_metrics = True
	#		 section_first_line = True
	#	 elif line.startswith('##Metrics global'):
	#		 remove_metrics = True
	#		 in_detection_metrics = False
	#		 section_first_line = True
	#	 elif line == "":
	#		 in_detection_metrics = False
	#		 in_remove_metrics = False
	#	 else:
	#		 if in_detection_metrics:
	#			 if section_first_line:
	#				 line_fields = line[1:].split("\t")[1:]
	#				 detection_categories = line_fields
	#				 section_first_line = False
	#			 else:
	#				 line_fields = line.split("\t")
	#				 detection_data.append({
	#						  'name': line_fields[0],
	#						  'data': map(int, line_fields[1:])
	#				 })
	#		 elif in_remove_metrics:
	#			 if section_first_line:
	#				 line_fields = line[1:].split("\t")
	#				 remove_categories = [category.lower().replace(" ", "_") for category in line_fields]
	#				 section_first_line = False
	#			 else:
	#				 for idx, val in enumerate(line.split("\t")):
	#					 remove_data[remove_categories[idx]] = int(val)
	# log_fh.close()

	# Write

	# FH_summary_tpl = open( os.path.join(CURRENT_DIR, "gene_placement.html") )
	# FH_summary_out = open( summary_file, "w" )
	# for line in FH_summary_tpl:
	#	 if "###DETECTION_CATEGORIES###" in line:
	#		 line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(detection_categories) )
	#	 elif "###DETECTION_DATA###" in line:
	#		 line = line.replace( "###DETECTION_DATA###", json.dumps(detection_data) )
	#	 elif "###LINE_NAME###" in line:
	#		 line = line.replace( "###LINE_NAME###", json.dumps(row) )
	#	 FH_summary_out.write( line )

	# FH_summary_out.close()
	# FH_summary_tpl.close()
##################### Faire un tableau de sortie : qui sera mon output: dico de dico ######################
#Function for create dico which contains all output_file
#dico: keys = titre colonne, value = contenue colonne 
def result_file(output):

	file  = open(output, "r")
	line = file.readline()

	global dico
	print("entete" +"_"+ line)
	tab_key = line.strip().split("\t")
	tab_col = []
	dico_2 = {}
	for i in tab_key:
		if i not in dico:
			dico[i] = []

			print("Voila le tableau")
		dico_2[i] = i
		#print(dico)

	for i in dico:
		if i in dico_2:
			tab_col.append(i)

#print(line.split(" "))
#print("test1")
	print("tab col : ", len(tab_col))
	for line in file:
		tab = line.strip().split("\t")
		print(line)
		e = 0
		for key in dico:
			print("e :", e)
			if e < len(tab_col) and  key == tab_col[e]:
				if e < len(tab):
					dico[tab_col[e]].append(tab[e])
				else:
					dico[tab_col[e]].append(":")
				e += 1
			else:
				dico[key].append(":")

	file.close()


def f2(output):
	file  = open(output, "r")
	line = file.readline()

	global dico
	dico[output] = {}
	print("entete" +"_"+ line)
	tab_key = line.strip().split("\t")
	for i in tab_key:
		if i not in dico[output]:
			dico[output][i] = []


	for line in file:
		tab = line.strip().split("\t")
		#print(line)
		e = 0
		for key in tab_key:
			dico[output][key].append(tab[e])
			e += 1

	file.close()
########### Dico fichier finale  ###################
#### Methode de Maria
def f3(output):
	global dict_tsv
	#dict_tsv = OrderedDict()

	FH = open(output,"rt")

	keys = FH.readline().strip().split()

	for line in FH:   
		for idx,value in enumerate(line.split()):
			if idx == 0:
				seq_id = value
				if seq_id not in dict_tsv:
					dict_tsv[seq_id] = OrderedDict()
			else :
				annot = keys[idx]
				dict_tsv[seq_id][annot] = value
	FH.close()

#########   Dico marker ###############
def f3marker(output):
	global dict_tsv1
	#dict_tsv = OrderedDict()

	FH1 = open(output,"rt")

	keys = FH1.readline().strip().split()

	for line in FH1:   
		for idx,value in enumerate(line.split()):
			if idx == 0:
				seq_id = value
				if seq_id not in dict_tsv1:
					dict_tsv1[seq_id] = OrderedDict()
			else :
				annot = keys[idx]
				dict_tsv1[seq_id][annot] = value
	FH1.close()
######### Dico reaction #################
def f3reaction(output):
	global dict_tsv1
	#dict_tsv = OrderedDict()

	FH1 = open(output,"rt")

	keys = FH1.readline().strip().split()

	for line in FH1:   
		for idx,value in enumerate(line.split()):
			if idx == 0:
				seq_id = value
				if seq_id not in dict_tsv1:
					dict_tsv1[seq_id] = OrderedDict()
			else :
				annot = keys[idx]
				dict_tsv1[seq_id][annot] = value
	FH1.close()

############ Fonction pour séparer les fichier marker et reaction ###########
def tsvParse(output, c, v):

	tsv_parse = "out_function"
	#tsv_parse = "Markers"+"_"+args.output
	# je crée un dataframe
	df = pd.read_csv(output, header=0, delimiter='\t')
	df_other = df[["sequence","metadata_NSTI", "16S_rRNA_Count"]]
	# Je suprime la partie marker_count
	al = df.drop(columns=["metadata_NSTI", "16S_rRNA_Count"])
	#print(df.drop(df.columns[[0, 1]], axis=1))
	print(al)
	# J écris dans un nouveau fichier
	al.to_csv(tsv_parse+str(v)+".tsv", "\t", index=None)
	df_other.to_csv("marker"+str(c)+".tsv", sep="\t", index=None)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
	#Categorie
	Categorie = ['16S', 'ITS', '18S']
	# Table to use for prediction
	TRAIT_OPTIONS = ['16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM', 'PHENO']
	 # Table to use for prediction
	FUNGI_OPTIONS = ['EC']
	# prediction method to use 
	HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

	#return value pair : 1er = indice (i) , 2em = value (v)
	find = False
	#c prend la chaine de caractére ajouter
	c = ""
	for i, v in enumerate(sys.argv):
		if sys.argv[i] == "-c" or sys.argv[i] == "--categorie":
			c = sys.argv[i+1]
			if sys.argv[i+1] == "16S":
				find = True	  

	# Manage parameters
	parser = argparse.ArgumentParser( description='predict gene family for OTU' )

	print("partie4")
	# Inputs
	#Input for categorie
	group_input = parser.add_argument_group( 'Inputs' )

	group_input.add_argument('-c', '--categorie',choices=Categorie, help='Specifies which categorie 16S or ITS, 18S')

	if find:
		group_input.add_argument('-i', '--in_trait',  help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')

	else:
		group_input.add_argument('-i', '--observed_trait_table', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')

	#Pour avoir plusieur parametres   
	# tab_arg = []
	# for i, v in enumerate(sys.argv):
	#	  if sys.argv[i] == "-i":
	#		 ind = 0
	#		 for e in range(i+1, len(sys.argv)):
	#			 if sys.argv[e][0] == "-" :
	#				 ind = e
	#				 print(i+1, ind, len(sys.argv), sys.argv[e])
	#				 break

	#		 if ind != 0:
	#			 tab_arg = sys.argv[i+1:ind]
	#			 print(tab_arg)
	#			 break

	#group_input.add_argument('-t', '--table', help='Specifies which table choices, ')

	#les parametres intermediares
	#group_input.add_argument('-i', '--in_trait', type=str.upper, choices=TRAIT_OPTIONS, help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')
	#group_input.add_argument('-i', '--in_trait',  help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')

	group_input.add_argument('-t', '--tree', metavar='PATH', required=True, type=str, help='The full reference tree in newick format containing both study sequences (i.e. ASVs or OTUs) and reference sequences.')

	group_input.add_argument('-o', '--output', metavar='PATH', type=str, required=True, help='Output table with predicted abundances per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')

	#group_input.add_argument('--observed', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')
	#group_input.add_argument('--observed_trait_table', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')


	group_input.add_argument('--chunk_size', default=500, type=int, help='Number of functions to run at a time on one processor. Note that you should consider how many ''processes you have specified before changing this ''option. E.g. if you specify the chunk_size to be ''the total number of functions, 1 processor will ''be used even if you specified more so the job will ''be substantially slower (default: %(default)d).')

	group_input.add_argument('-m', '--hsp_method', default='mp',choices=HSP_METHODS, help='HSP method to use.' +'"mp": predict discrete traits using max parsimony. ''"emp_prob": predict discrete traits based on empirical ''state probabilities across tips. "subtree_average": ''predict continuous traits using subtree averaging. ' '"pic": predict continuous traits with phylogentic ' 'independent contrast. "scp": reconstruct continuous ''traits using squared-change parsimony (default: ''%(default)s).')

	group_input.add_argument('-p', '--processes', default=1, type=int, help='Number of processes to run in parallel (default: ' '%(default)d).')

	group_input.add_argument('--verbose', default=False, action='store_true',help='If specified, print out wrapped commands and other ''details to screen.')

	group_input.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)

	#group_input.add_argument('-b', '--biom_file', metavar='PATH', required=True, type=str, help='Biom file.')

  
	# output
	group_output = parser.add_argument_group( 'Outputs' )

	#group_output.add_argument('-s','--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )	

	group_output.add_argument('-n', '--calculate', default=False, action='store_true', help='Calculate NSTI and add to output file.')

	group_output.add_argument('--check', default=False, action='store_true', help='Check input trait table before HSP.')

	group_output.add_argument('--seed', default=100, type=int, help='Seed to make output reproducible, which is ''necessary for the emp_prob method ' '(default: %(default)d).')

	group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

	print("partie5")
	args = parser.parse_args()
	"""
	@tab_arg: table of value wich contain trait-table (EC,KO,PFAM...)
	@ i : indice of value
	@ v : value
	@ note : values must be separated by commas ","
	"""
	tab_arg = []
	print(tab_arg)
	if find: 
		tab_arg = args.in_trait.split(",")
		tab_arg.append("16S")
	else :
		tab_arg = args.observed_trait_table.split(",")
		#tab_arg.append(c)
		#tab_arg.append("ITS")
		#tab_arg.append("18S")
	print(tab_arg)
	#Ajouter 16S directement 

	## Verifier le contenu de tab_arg avant de lancer ##
	print(tab_arg)
	for i in range(len(tab_arg)):
		print("lllllllll")
		print(tab_arg[i] not in TRAIT_OPTIONS)
		print("oooooooooo")
		print(tab_arg[i])
		if find == True and  tab_arg[i] not in TRAIT_OPTIONS:
			print(tab_arg[i])
			print("You choice in " + str(TRAIT_OPTIONS))
			exit() 
		if find == False and tab_arg[i] not in FUNGI_OPTIONS:
			print("You choice in " + str(FUNGI_OPTIONS))
			exit()

	####Declare error output ###
	stderr = "hsp.stderr"
	### Declare hsp methods ###
	hsp16SMet = ""
	hspITSMet = ""
	print("partie6")
	# Process 
	#tab_marker = []
	tab_reaction = []
	tab_files = []
	try:	 
		Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
		# Commands execution
		print("partie7", tab_arg)
		#Test = False
		keep = ""
		for i, v in enumerate(tab_arg+[c]):
			##### Debut bon script ####	  
			# if args.categorie == "16S":
			#	 #hsp16_cmd = hsp16S(hsp16SMet ,"16S", v, args.tree, "16S_count_table_"+str(i)+"_"+args.output, args.calculate, stderr).submit(args.log_file)
			#	 hsp16_cmd = hsp16S(hsp16SMet ,"16S", v, args.tree, "reaction_"+str(i)+"_"+args.output, args.calculate, stderr).submit(args.log_file)

			# #Separer les deux fichier (qu on les distingue)
			# elif args.categorie == "ITS" and v == "EC": 
			#	  hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "reaction_count_table"+":"+args.output, EC, args.calculate, stderr).submit(args.log_file) 
			# elif args.categorie == "ITS" : 
			#	 hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "ITS_count_table"+":"+args.output, ITS, args.calculate, stderr).submit(args.log_file)  
			# elif args.categorie == "18" and v == "EC":
			#	 hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "reaction_count_table_"+str(i)+"_"+args.output, EC_18, args.calculate, stderr).submit(args.log_file)		  
			# elif args.categorie == "18S":
			#	 hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "18S_count_table_"+str(i)+"_"+args.output, fungi_18S, args.calculate, stderr).submit(args.log_file)  
			# #Sortir les 16S ici #tab arg prend pas la sortie 16S (modifier l'agr)
			# # tab_marker.append("ITS_count_table"+":"+args.output)
			# # tab_reaction.append("ITS_count_table"+":"+args.output)
			# tab_files.append("16S_count_table_"+str(i)+"_"+args.output) 
			# #tab_files.append("fichier_"+str(i)+"_"+args.output)  
			####### Fin bon script ############
			if args.categorie == "16S":
				hsp16_cmd = hsp16S(hsp16SMet ,"16S", v, args.tree, "fichier_"+str(i)+"_"+args.output, args.calculate, stderr).submit(args.log_file)
				
			elif args.categorie == "ITS" and v == "EC": 
				hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, EC, args.calculate, stderr).submit(args.log_file) 

			elif args.categorie == "ITS" : 
				hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, ITS, args.calculate, stderr).submit(args.log_file)  
			elif args.categorie == "18" and v == "EC":
				hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, EC_18, args.calculate, stderr).submit(args.log_file)		  
			elif args.categorie == "18S":
				hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, fungi_18S, args.calculate, stderr).submit(args.log_file)  
			keep += "_" + str(v)
			tab_files.append("fichier_"+str(i)+"_"+args.output) 
		############ Pour marker###########
		# dict_tsv1 = {}
		
		# for output in tab_marker:
		#	 f4(output)


		# fout1 = open("Masorti"+args.output, "w")
		# print("test1 Markeur")
		# # for k, v in dict_tsv.items():
		# #	 fout.write(str(k) + "\t" + str(v) + "\n" )
		# # fout.close()
		# header ="sequence"
		# print(dict_tsv1, "dict_tsv-------------------")
		# for seq_id in dict_tsv1:
		#	 print("aaaaaa#1", seq_id)
		#	 if header == "sequence":
		#		 header = header + "\t" + "\t".join(dict_tsv1[seq_id].keys()) + "\n"
		#		 fout1.write(header)
		#		 print("aaaaaa#lllll", seq_id)
		#	 line=seq_id
		#	 for annot in dict_tsv1[seq_id]:
		#		 line=line + "\t" + dict_tsv1[seq_id][annot]
		#	 fout1.write(line + "\n")  
		#	 print("45545#1", seq_id)
		# fout1.close()
		# exit()
		########## Fin marker ###################
		
		# #write_summary(args.html, args.output)

		# ### Methode Maria:
		dict_tsv ={}
		for output in tab_files:
			f3(output)

		fout = open(args.output, "w")
		print("test1")
		header ="sequence"
		print(dict_tsv, "dict_tsv-------------------")
		for seq_id in dict_tsv:
			print("aaaaaa#1", seq_id)
			if header == "sequence":
				header = header + "\t" + "\t".join(dict_tsv[seq_id].keys()) + "\n"
				fout.write(header)
				print("aaaaaa#lllll", seq_id)
			line=seq_id
			for annot in dict_tsv[seq_id]:
				line=line + "\t" + dict_tsv[seq_id][annot]
			fout.write(line + "\n")  
			print("45545#1", seq_id)
		fout.close()


		# Apelle à la fonction de parse du tsv
		tsvParse(args.output, c, keep)
		#tsvParse("fichier_"+str(i)+"_"+args.output)
		print("tsvParse FINI")

		######################################
		# dict_tsv = {}
		
		# for output in tab_marker:
		#	 f4(output)

		# #print(dico)
		# #print(result_file)
		# print(tab_marker)

		# fout2 = open(args.output, "w")
		# print("test2")
		# # for k, v in dict_tsv.items():
		# #	 fout.write(str(k) + "\t" + str(v) + "\n" )
		# # fout.close()
		# header ="sequence"
		# print(dict_tsv, "dict_tsv-------------------")
		# for seq_id in dict_tsv:
		#	 print("aaaaaa#1", seq_id)
		#	 if header == "sequence":
		#		 header = header + "\t" + "\t".join(dict_tsv[seq_id].keys()) + "\n"
		#		 fout2.write(header)
		#		 print("aaaaaa#lllll", seq_id)
		#	 line=seq_id
		#	 for annot in dict_tsv[seq_id]:
		#		 line=line + "\t" + dict_tsv[seq_id][annot]
		#	 fout2.write(line + "\n")  
		#	 print("45545#1", seq_id)
		# fout2.close()
		# exit()
			#else :
					#hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, v, args.calculate, stderr).submit(args.log_file)		 
			
	
		#hsp_16S(hsp_16_Met, args.in_trait, args.tree, args.output, args.observed_trait_table, args.calculate, stderr)
		#hsp16S(hsp16SMet, args.in_trait, args.tree,args.output,args.calculate, stderr).submit(args.log_file)
		#
		#print("partie8")
		#hsp_EC(hsp_EC_Met, args.tree, args.output, args.observed_trait_table, args.calculate, stderr)

		
		#hspITS(hspITSMet, args.tree, args.output,  args.observed_trait_table, args.calculate, stderr).submit(args.log_file)
		#hspITS(hspITSMet, args.tree, args.output,  args.observed, args.calculate, stderr)

		#print("partie9")	
		#write_summary(args.html, "output", args.biom_file)

		# Result_out = {args.output}
		# print(Result_out)
		# Result_out

		################
		# dict_tsv = {}
		
		# for output in tab_files:
		#	 f3(output)

		# # #print(dico)
		# # #print(result_file)
		# # print(tab_files)

		# # # fout = open("Mous.txt", "w")

		# # # i = 0
		# # # p = list(dico[tab_files[0]].keys())[0]
		# # # size = len(dico[tab_files[0]][p])


		# # # for fi in dico:
		# # #	 fout.write("-:" + fi + ":-\n")
		# # #	 fout.write( "\t".join(list(dico[fi].keys())) + "\n" )

		# # #	 #print("\t".join(list(dico[fi].keys())))
		# # #	 row = 0
		# # #	 size = len(dico[fi][dico[fi].keys()[0]])
		# # #	 while row < size :
		# # #		 for col in dico[fi]:
		# # #			 if row < len(dico[fi][col]):
		# # #				 #print(dico[fi][col])
		# # #				 fout.write(str(dico[fi][col][row]) + "\t")

		# # #		 fout.write("\n")

		# # #		 row += 1
		# # # fout.write("\n")

		# # # fout.close()
		# # # exit()
		# # # ### Test dico
		# # # ft = open("tttt.txt", "w")
		# # # ft.write(str(dico))
		# # # ft.close()

		
	finally:
		print("Partie Finale ")