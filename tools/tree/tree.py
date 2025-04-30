#!/usr/bin/env python3

__author__ = 'Ta Thi Ngan - SIGENAE/GABI & Maria Bernard - SIGENAE/GABI'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse
import json
from numpy import median

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']
os.environ['MAFFT_BINARIES'] = ""

from frogsUtils import *
from frogsSequenceIO import *
from frogsBiom import BiomIO
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Mafft(Cmd):
    """
    @summary: Mafft denovo multiple alignment of observed candidate sequences.
    @see: http://mafft.cbrc.jp/alignment/software/
    """
    def __init__(self, mafftMet, fasta, aligned, thread, stderr):
        """
        @param mafftMet: [str] Mafft method option.
        @param fasta: [str] Path to input fasta file.
        @param aligned: [str] Path to store resulting alignment file.
        @param stderr: [str] Path to temporary Mafft stderr output file
        @param thread: [int] number of cpu to use.
        """
        Cmd.__init__(self,
                 'mafft',
                 'Mafft multiple alignment.',
                  "" + mafftMet + "--thread "+ str(thread) +" "+ fasta +" > "+aligned +' 2> ' + stderr,
                  "--version")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').split()[0].strip()

class FastTree(Cmd):
    """
    @summary: reconstruction of phylogenetic tree with the GTR model for nucleotic evolution.
    @see: http://www.microbesonline.org/fasttree/
    """
    def __init__(self, align, output, fasttree_stderr):
        """
        @param align: [str] Path to input alignment file.
        @param output: [str] path to store resulting tree file.
        @param fasttre_stderr: [str] Path to temporary FastTree stderr output file
        """
        Cmd.__init__( self,
                      'FastTree',
                      'reconstruction a phylogenetic tree',
                      "-nt -gtr " + align+" > "+ output+ ' 2> ' + fasttree_stderr,
                      "-expert")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').split()[4].strip()


class RootTree(Cmd):
    """
    @summary: root tree with phangornm midpoint
    @see: https://cran.r-project.org/web/packages/phangorn/phangorn.pdf
    """
    def __init__(self, in_tree, out_tree):
        """
        @param in_tree: [str] Path to input tree file (Newick format).
        @param out_tree: [str] path to output rooted tree file (Newick format)
        """
        Cmd.__init__( self,
                    "root_tree.R",
                    "root newick tree with phangorn R package midpoint function.",
                    in_tree + " " + out_tree,
                    "-v")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def get_fasta_nb_seq( fasta_file ):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file to process.
    @return: [int] The number of sequences.
    """
    return sum(1 for _ in (FastaIO(fasta_file)))

def get_methods_mafft(seqs):
    """
    @summary: return mafft options
    @param seqs : [str] Path to input fasta file.
    @return [str] :  "--maxiterate 1000 --globalpair " if number of sequences is less than or equal to 200
                     "--maxiterate 1000 " if number of sequences is greater than 200 and less than 500
                     "--auto " if number of sequences is greater than 500
    """
    n=sum(1 for _ in (FastaIO(seqs)))
    if n<=200:
        return "--maxiterate 1000 --globalpair "
    elif 500>n>200:
        return "--maxiterate 1000 "
    else:
        return "--auto "
        
def write_summary( summary_file, fasta_in, align_out, biomfile, treefile):
    """
    @summary: Writes the process summary in one html file.
    @param summary_file: [str] path to the output html file.
    @param align_out: [str] path to the fasta file of unaligned ASV
    @param biomfile: [str] path to the input BIOM file.
    @param treefile: [str] path to the Newick file.
    """
    # to summary ASVs number && abundances number               
    summary_info = {
       'asv_kept' : 0,
       'asv_removed' : 0,
       'abundance_kept' : 0,
       'abundance_removed' : 0       
    }
    number_asv_all = 0
    number_abundance_all = 0
    # to detail removed ASV
    removed_details_categories =["Taxonomic Information", "Abundance Number", "% with abundance total", "Sequence length"]
    removed_details_data =[]
    
    # to build one metadata for tree view
    dic_asv={}
    list_asv_all=list()
    list_out_tree=[]

    biom=BiomIO.from_json(biomfile)
    treefile = open(treefile, "rt")
    newick = treefile.read().strip()

    # record nb ASV and abundance
    for asv in FastaIO(fasta_in):
        list_asv_all.append(asv.id)
        number_asv_all +=1
        number_abundance_all += biom.get_observation_count(asv.id)

    # record details about removed ASV
    if align_out is not None:
        for asv in FastaIO(align_out):
            summary_info['asv_removed'] +=1
            summary_info['abundance_removed'] += biom.get_observation_count(asv.id)
            
            # to built one table of ASVs out of phylogenetic tree
            taxonomy=""
            if biom.has_metadata("taxonomy"):
                taxonomy = ";".join(biom.get_observation_metadata(asv.id)["taxonomy"]) if issubclass(biom.get_observation_metadata(asv.id)["taxonomy"].__class__,list) else str(biom.get_observation_metadata(asv.id)["taxonomy"])
            elif biom.has_metadata("blast_taxonomy"): 
                taxonomy = ";".join(biom.get_observation_metadata(asv.id)["blast_taxonomy"]) if issubclass(biom.get_observation_metadata(asv.id)["blast_taxonomy"].__class__,list) else str(biom.get_observation_metadata(asv.id)["blast_taxonomy"])
            abundance=biom.get_observation_count(asv.id)
            percent_abundance=abundance*100/(float(number_abundance_all))
            length=len(asv.string)
            info={"name": asv.id, "data": [taxonomy, abundance, percent_abundance, length]}
            removed_details_data.append(info)
            list_out_tree.append(asv.id)

    # improve tree view by adding taxonomy information
    list_in_tree=[item for item in list_asv_all if item not in list_out_tree]  
    for asv in list_in_tree:
        tax=None
        if biom.has_metadata("taxonomy"):
            tax=" ".join(biom.get_observation_metadata(asv)["taxonomy"]) if issubclass(biom.get_observation_metadata(asv)["taxonomy"].__class__, list) else str(biom.get_observation_metadata(asv)["taxonomy"])
        elif biom.has_metadata("blast_taxonomy"):
            tax=" ".join(biom.get_observation_metadata(asv)["blast_taxonomy"]) if issubclass(biom.get_observation_metadata(asv)["blast_taxonomy"].__class__, list) else str(biom.get_observation_metadata(asv)["blast_taxonomy"])
        if tax :
            newick=newick.replace(asv + ":", asv + " " + tax + ":")
    
    # finalize summary
    summary_info['asv_kept'] = number_asv_all - summary_info['asv_removed']
    summary_info['abundance_kept'] = number_abundance_all - summary_info['abundance_removed']
    
    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "tree_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    for line in FH_summary_tpl:
        if "###HEIGHT###" in line:
            line = line.replace( "###HEIGHT###", json.dumps(summary_info['asv_kept']*11+166))
        if "###NEWICK###" in line:
            line = line.replace( "###NEWICK###", newick)
        if "##REMOVED_DETAILS_CATEGORIES###" in line:
            line = line.replace( "###REMOVED_DETAILS_CATEGORIES###", json.dumps(removed_details_categories) )
        elif "###REMOVED_DETAILS_DATA###" in line:
            line = line.replace( "###REMOVED_DETAILS_DATA###", json.dumps(removed_details_data) )
        elif "###SUMMARY###" in line:
            line = line.replace( "###SUMMARY###", json.dumps(summary_info) )
        elif '<div id="ASVs-fail" style="display:none;">' in line:
            if summary_info['asv_removed']!=0:
                line = line.replace( 'style="display:none;"', '' )
        elif "###FROGS_VERSION###" in line:
            line = line.replace( "###FROGS_VERSION###", "\""+str(__version__)+"\"" )
        elif "###FROGS_TOOL###" in line:
            line = line.replace( "###FROGS_TOOL###", "\""+ os.path.basename(__file__)+"\"" )
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
    parser.add_argument('--version', action='version', version=__version__ )
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument('--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )
    
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--input-fasta', required=True, help='Path to input FASTA file of ASV seed sequences. Warning: FROGS Tree is only working on less than 10000 sequences!' )
    group_input.add_argument('--input-biom', required=True, help='Path to the abundance BIOM file.' )
        
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--output-tree', default='tree.nwk', help="Path to store resulting Newick tree file. (format: nwk) [Default: %(default)s]" )
    group_output.add_argument('--html', default='tree.html', help="The HTML file containing the graphs. [Default: %(default)s]" )    
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)
    
    ### Temporary files
    tmpFiles=TmpFiles(os.path.split(args.output_tree)[0])
    filename_prefix = ".".join(os.path.split(args.input_fasta)[1].split('.')[:-1])
    
    # alignment temporary files
    stderr = tmpFiles.add("mafft.stderr")
    align= tmpFiles.add('mafft_aligned.fasta')   
    # if we want to add alignment method that do not keep necessarily all ASV, such as pynast when supported in FROGS
    align_out=None 
    
    # fastree temporary files
    fasttree=tmpFiles.add("fasttree.nwk")
    fasttree_stderr=tmpFiles.add("fasttree.stderr")
    
    # Process 
    try:        
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        nb_seq = get_fasta_nb_seq(args.input_fasta)
        biom = BiomIO.from_json(args.input_biom)
        if nb_seq > len(biom.rows):
            raise_exception( Exception("\n\n#ERROR : Your fasta input file contains more ASV than your biom file.\n\n"))
        Logger.static_write(args.log_file, "Number of input ASVs sequences: " + str(nb_seq) + "\n\n")
        if nb_seq >10000:
            raise_exception( Exception( "\n\n#ERROR : FROGS Tree is only working on less than 10 000 sequences!\n\n" ))
        
        # alignment step
        mafftMet=get_methods_mafft(args.input_fasta)
        Mafft(mafftMet, args.input_fasta, align, args.nb_cpus, stderr).submit( args.log_file )

        # tree contruction step
        FastTree(align, fasttree, fasttree_stderr).submit( args.log_file )

        # rooting tree step
        RootTree(fasttree, args.output_tree).submit(args.log_file)

        # summarize resultats in HTML output 
        write_summary( args.html, args.input_fasta, align_out, args.input_biom, args.output_tree)
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
