#!/usr/bin/env python3
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

import os
import sys
import argparse
import json
import re
from numpy import median

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR #
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsSequenceIO import *
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES 
#
##################################################################################################################################################
class picrust2_place_seqs(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: place_seqs.py program placed Fasta sequence on the tree Arbre .
    @see: https://github.com/picrust/picrust2/wiki
    """
    def __init__(self, study_fasta, out_tree,min_align,ref_dir, stdout):
        """
        @param mafftMet: [str] picrust2 method option.
        @param fasta: [str] Path to input fasta file.
        @param out_tree: [str] Path to store resulting tree file.
        @param stderr: [str] Path to temporary picrust2 stderr output file
        @param thread: [int] number of cpu to use.
        """
        Cmd.__init__(self,
                 'place_seqs.py',
                 'place OTUs on tree.',
                  "--study_fasta "+ str(study_fasta) +" --out_tree "+ str(out_tree) +" --min_align "+ str(min_align) +" --ref_dir "+ str(ref_dir) +' 2> ' + stdout,
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
def restricted_float(in_arg):
    """
    @summary: Custom argparse type to force an input float to be between 0 and 1.
    """
    try:
        in_arg = float(in_arg)
    except ValueError:
        raise argparse.ArgumentTypeError(in_arg + " is not a floating-point "
                                         "literal (i.e. not a proportion)")

    if in_arg < 0.0 or in_arg > 1.0:
        raise argparse.ArgumentTypeError(in_arg + "is not in range 0.0 - 1.0")
    return in_arg

def get_fasta_nb_seq(fasta_file):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file to process.
    @return: [int] The number of sequences.
    @ Adapt fasta file to picrust2 input
    """
    return sum(1 for _ in (FastaIO(fasta_file)))

def convert_fasta(fasta_file):
    """
    @summary: Change fasta headers to be compatible with picrust2
    """

    FH_input = FastaIO(fasta_file)
    FH_output = FastaIO( "sout.fasta","wt" )
    for record in FH_input:
        record.id = record.id
        record.description = None
        FH_output.write(record)
    FH_output.close()
   
def excluded_sequence(file_tree, fasta_file, out_file):
    """
    @summary: Returns the excluded sequence.
    @param fasta_file: [str] Path to the fasta file to process.
    @param tree_file: [str] Path to the tree file to process.
    @return: [int] The file of no aligned sequence.
    """
    file = open(file_tree, "r")
    line = file.readline()
    list_cluster = re.findall("(Cluster_[0-9]+)", line)
    file.close()

    FH_input = FastaIO(fasta_file)
    FH_output = FastaIO(out_file, "wt")

    for record in FH_input:
        if record.id not in list_cluster:
            FH_output.write(record)
    FH_input.close()
    FH_output.close()
 
def test(file_tree, biom):
    file = open(file_tree, "r")
    line = file.readline()
    #List of cluster
    list_cluster = []
    #Boucle sur le fichier tree
    #Je splite sur (,) , regex sur le cluster et récupérer le groupe1
    #Je parcour ligne par ligne
    while line: 
        for i, v in enumerate(line.split(",")):
            group = re.search("(Cluster_[0-9]+)", v)
            if group:
                ide = group.group(1)
                list_cluster.append(ide)

        line = file.readline() 

    file.close()

    file_biom = open(biom, "r")

    file_out1 = open(biom.split(".")[0] + "_removed.tsv", "w")
    file_out2 = open(biom.split(".")[0] + "_kept.tsv", "w")

    line = file_biom.readline()
    dico_biom = json.loads(line.strip())
    print(dico_biom)

    for i, value in enumerate(dico_biom["rows"]) :
        if dico_biom["rows"][i]["id"] in list_cluster :
            file_out1.write(dico_biom["rows"][i]["id"]+"\n")
        else:
            file_out2.write(dico_biom["rows"][i]["id"]+"\n")


    file_biom.close()
    file_out1.close()
    file_out2.close()


def write_summary( summary_file, fasta_in, align_out, biomfile, treefile ):
	
    """
    @summary: Writes the process summary in one html file.
    @param summary_file: [str] path to the output html file.
    @param align_out: [str] path to the fasta file of unaligned OTU
    @param biomfile: [str] path to the input BIOM file.
    @param treefile: [str] path to the Newick file.
    """
    # to summary OTUs number && abundances number               
    summary_info = {
       'otu_kept' : 0,
       'otu_removed' : 0,
       'abundance_kept' : 0,
       'abundance_removed' : 0       
    }
    number_otu_all = 0
    number_abundance_all = 0
    # to detail removed OTU
    removed_details_categories =["Taxonomic Information", "Abundance Number", "% with abundance total", "Sequence length"]
    removed_details_data =[]
    
    # to build one metadata for tree view
    dic_otu={}
    list_otu_all=list()
    list_out_tree=[]

    biom=BiomIO.from_json(biomfile)
    treefile = open(treefile, "r")
    newick = treefile.read().strip()

    # record nb OTU and abundance
    for otu in FastaIO(fasta_in):
        list_otu_all.append(otu.id)
        number_otu_all +=1
        number_abundance_all += biom.get_observation_count(otu.id)

    # record details about removed OTU
    if align_out is not None:
        for otu in FastaIO(align_out):
            summary_info['otu_removed'] +=1
            summary_info['abundance_removed'] += biom.get_observation_count(otu.id)
            
            # to built one table of OTUs out of phylogenetic tree
            taxonomy=""
            if biom.has_metadata("taxonomy"):
                taxonomy = ";".join(biom.get_observation_metadata(otu.id)["taxonomy"]) if issubclass(biom.get_observation_metadata(otu.id)["taxonomy"].__class__,list) else str(biom.get_observation_metadata(otu.id)["taxonomy"])
            elif biom.has_metadata("blast_taxonomy"): 
                taxonomy = ";".join(biom.get_observation_metadata(otu.id)["blast_taxonomy"]) if issubclass(biom.get_observation_metadata(otu.id)["blast_taxonomy"].__class__,list) else str(biom.get_observation_metadata(otu.id)["blast_taxonomy"])
            abundance=biom.get_observation_count(otu.id)
            percent_abundance=abundance*100/(float(number_abundance_all))
            length=len(otu.string)
            info={"name": otu.id, "data": [taxonomy, abundance, percent_abundance, length]}
            removed_details_data.append(info)
            list_out_tree.append(otu.id)

    # improve tree view by adding taxonomy information
    list_in_tree=[item for item in list_otu_all if item not in list_out_tree]  
    for otu in list_in_tree:
        tax=None
        if biom.has_metadata("taxonomy"):
            tax=" ".join(biom.get_observation_metadata(otu)["taxonomy"]) if issubclass(biom.get_observation_metadata(otu)["taxonomy"].__class__, list) else str(biom.get_observation_metadata(otu)["taxonomy"])
        elif biom.has_metadata("blast_taxonomy"):
            tax=" ".join(biom.get_observation_metadata(otu)["blast_taxonomy"]) if issubclass(biom.get_observation_metadata(otu)["blast_taxonomy"].__class__, list) else str(biom.get_observation_metadata(otu)["blast_taxonomy"])
        if tax :
            newick=newick.replace(otu + ":", otu + " " + tax + ":")
    
    # finalize summary
    summary_info['otu_kept'] = number_otu_all - summary_info['otu_removed']
    summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']
    
    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "FPStep1.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###HEIGHT###" in line:
            line = line.replace( "###HEIGHT###", json.dumps(summary_info['otu_kept']*11+166))
        if "###NEWICK###" in line:
            line = line.replace( "###NEWICK###", newick)
        if "##REMOVED_DETAILS_CATEGORIES###" in line:
            line = line.replace( "###REMOVED_DETAILS_CATEGORIES###", json.dumps(removed_details_categories) )
        elif "###REMOVED_DETAILS_DATA###" in line:
            line = line.replace( "###REMOVED_DETAILS_DATA###", json.dumps(removed_details_data) )
        elif "###SUMMARY###" in line:
            line = line.replace( "###SUMMARY###", json.dumps(summary_info) )
        elif '<div id="OTUs-fail" style="display:none;">' in line:
            if summary_info['otu_removed']!=0:
                line = line.replace( 'style="display:none;"', '' )
        FH_summary_out.write(line)
    FH_summary_out.close()
    FH_summary_tpl.close()
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":

    # Manage parameters
    parser = argparse.ArgumentParser( description='Phylogenetic tree reconstruction' )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )

    group_input.add_argument('-s', '--study_fasta', metavar='PATH', required=True,type=str, help='FASTA of unaligned study sequences.')

    group_input.add_argument('-b', '--biom_file', metavar='PATH', required=True, type=str, help='Biom file.')

    group_input.add_argument('-p', '--processes', type=int, default=1, help='Number of processes to run in parallel (default: ''%(default)d). Note that this refers to ''multithreading rather than multiprocessing when ''running EPA-ng and GAPPA.')

    group_input.add_argument('-r', '--ref_dir', metavar='PATH', type=str, required=True, help='Directory containing reference sequence files ''(default: %(default)s). Please see the online ''documentation for how to name the files in this ''directory in order to use custom reference files.')

    group_input.add_argument('--min_align', type=restricted_float, default=0.8, help='Proportion of the total length of an input query ''sequence that must align with reference sequences. ''Any sequences with lengths below this value after ''making an alignment with reference sequences will ''be excluded from the placement and all subsequent ''steps. (default: %(default)d).')

    group_input.add_argument('--chunk_size', type=int, default=5000, help='Number of query seqs to read in at once for EPA-ng ''(default: %(default)d).')

    group_input.add_argument('--verbose', default=False, action='store_true', help='If specified, print out wrapped commands and other ''details to screen.')

    group_input.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)

  
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-m','--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )    

    group_output.add_argument('-o', '--out_tree', metavar='PATH', required=True, type=str, help='Name of final output tree.')

    group_output.add_argument('--intermediate', metavar='PATH', type=str, default=None, help='Output folder for intermediate files (will be ''deleted otherwise).')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

    args = parser.parse_args()

    stderr = "picrust2.stderr"

    # Process 
    try:     
        print("\n\n\n--------------")
        os.system("pwd")
        #place_seqs.py --study_fasta --min_align  --out_tree --ref_dir --threads
        convert_fasta(args.study_fasta)
        picrust2_place_seqs("sout.fasta", args.out_tree, args.min_align, args.ref_dir, stderr).submit( args.log_file )
        #picrust2_place_seqs("sout.fasta", args.out_tree, args.min_align, stderr)
        print("Partie picrust fini")
        excluded_sequence(args.out_tree, "sout.fasta", "excluded.tsv")
        #PICRUSt2("sout.fasta", args.out_tree, args.min_align, args.ref_dir, stderr).submit( args.log_file )
        print("Partie 2 ")
        

        #RootTree(fasttree, args.out_tree).submit(args.log_file)

       # write_summary("sorti.html", "sout.fasta", "soutuiui.fasta", args.biom_file, args.out_tree )

        write_summary( args.html, "sout.fasta", "excluded.fasta", args.biom_file, args.out_tree)
    finally:
        print("Partie Finale ")






        
