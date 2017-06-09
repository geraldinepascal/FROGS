#!/usr/bin/env python2.7
#
# Copyright (C) 2017 INRA
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

__author__ = ' Ta Thi Ngan & Maria Bernard INRA - SIGENAE '
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs@toulouse.inra.fr'
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
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsUtils import *
from frogsSequenceIO import *
from frogsBiom import BiomIO
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Pynast(Cmd):
    """
    @summary: Multiple alignment of observed candidate sequences with aligned well-known sequences.
    @see: https://github.com/biocore/pynast
    """
    def __init__(self, in_candi_pynast, in_templ_pynast, min_len, pynast_align, pynast_fail, pynast_log):
        """
        @param in_candi_pynast: [str] Path to input fasta file.
        @param in_templ_pynast: [str] Path to template alignment file.
        @param min_len: [int] Minimum sequence length to include in NAST alignment.
        @param align: [str] Path to Pynast aligned sequences output.
        @param pynast_fail: [str] Path to Pynast failed aligned sequences output.
        @param pynast_log: [str] Path to Pynast log file.
        """
        Cmd.__init__( self,
                      'pynast',
                      'Pynast multiple alignment.',
                      "--input_fp " + in_candi_pynast + " --template_fp " + in_templ_pynast + " --min_len " + min_len + " --fasta_out "+ pynast_align + " --failure_fp " + pynast_fail + " --log_fp " + pynast_log + ' 2> /dev/null',
                      "--version")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').split()[2].strip()

class Mafft(Cmd):
    """
    @summary: Mafft denovo multiple alignment of observed candidate sequences.
    @see: http://mafft.cbrc.jp/alignment/software/
    """
    def __init__(self, mafftMet, fasta, aligned, thread):
        """
        @param mafftMet: [str] Mafft method option.
        @param fasta: [str] Path to input fasta file.
        @param aligned: [str] Path to store resulting alignment file.
        @param thread: [int] number of cpu to use.
        """
        Cmd.__init__(self,
                 'mafft',
                 'Mafft multiple alignment.',
                  "" + mafftMet + "--thread "+ str(thread) +" "+ fasta +" > "+aligned +' 2> /dev/null',
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
    def __init__(self, align, output):
        """
        @param align: [str] Path to input alignment file.
        @param output: [str] path to store resulting tree file.
        """
        Cmd.__init__( self,
                      'FastTree',
                      'reconstruction a phylogenetic tree',
                      "-nt -gtr " + align+" > "+ output+ ' 2> /dev/null',
                      "-expert")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').split()[4].strip()

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

def compute_min_sequence_length(seqs):
    """
    @summary:  return the minimum sequence length.
    @param seqs : [str] Path to input fasta file.
    """
    FH_in=FastaIO(seqs)
    min_length = min([len(record.string) for record in FH_in ])
    return str(min_length)

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
        
def write_summary( summary_file, pynast_fail, biomfile, treefile):
    """
    @summary: Writes the process summary in one html file.
    @param summary_file: [str] path to the output html file.
    @param pynast_fail: [str] path to the pynast_fail file.
    @param biomfile: [str] path to the input BIOM file.
    @param treefile: [str] path to the Newick file.
    """
    # to summary OTUs number && abundances number				
    detection_categories =["Taxonomic Information", "Abundance Number", "% with abundance total", "Sequence length"]
    table_otu_out =[]
    summary_info = {
       'number_otu_all' : 0,
       'otu_kept' : 0,
       'otu_removed' : 0,
       'number_abundance_all' : 0,
       'abundance_kept' : 0,
       'abundance_removed' : 0       
    }
    biom=BiomIO.from_json(biomfile)
    # to build one metadata for tree view
    dic_otu={}
    treefile = open(treefile, "r")
    newick = treefile.read().strip()
    list_otu_all=biom.get_observations_names()
    list_out_tree=[]
    for observation_name in biom.get_observations_names():
        summary_info['number_otu_all'] +=1
        summary_info['number_abundance_all'] += biom.get_observation_count(observation_name)

    if pynast_fail is not None:
        for otu in FastaIO(pynast_fail):
            summary_info['otu_removed'] +=1
            summary_info['abundance_removed'] += biom.get_observation_count(otu.id)
            
            # to built one table of OTUs out of phylogenetic tree
            taxonomy=""
            if biom.has_metadata("taxonomy") or biom.has_metadata("blast_taxonomy"):
                taxonomy=";".join(biom.get_observation_taxonomy( otu.id, "taxonomy" )) if biom.has_observation_metadata( 'taxonomy' ) else ";".join(biom.get_observation_taxonomy( otu.id, "blast_taxonomy" ))
            abundance=biom.get_observation_count(otu.id)
            percent_abundance=abundance*100/(float(summary_info['number_abundance_all']))
            length=len(otu.string)
            info={"name": otu.id, "data": [taxonomy, abundance, percent_abundance, length]}
            table_otu_out.append(info)
            list_out_tree.append(otu.id)

    list_in_tree=[item for item in list_otu_all if item not in list_out_tree]
    
    for otu in list_in_tree:
        tax=None
        if biom.has_metadata("taxonomy") or biom.has_metadata("blast_taxonomy"):
            tax=" ".join(biom.get_observation_taxonomy(otu, "taxonomy" )) if biom.has_observation_metadata( 'taxonomy' ) else " ".join(biom.get_observation_taxonomy(otu, "blast_taxonomy"))
        if tax :
            newick=newick.replace(otu + ":", otu + " " + tax + ":")
            
    summary_info['otu_kept'] = summary_info['number_otu_all'] - summary_info['otu_removed']
    summary_info['abundance_kept'] = summary_info['number_abundance_all'] - summary_info['abundance_removed']
    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "tree_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###HEIGHT###" in line:
            line = line.replace( "###HEIGHT###", json.dumps(summary_info['otu_kept']*11+166))
        if "###NEWICK###" in line:
            line = line.replace( "###NEWICK###", newick)
        if "###DETECTION_CATEGORIES###" in line:
            line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(detection_categories) )
        elif "###DETECTION_DATA###" in line:
            line = line.replace( "###DETECTION_DATA###", json.dumps(table_otu_out) )
        elif "###REMOVE_DATA###" in line:
            line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
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
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]" )

    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-i', '--input-otu', required=True, help='Path to input fasta file of OTU. Warning: FROGS Tree is only working on less than 10000 sequences!' )
    group_input.add_argument( '-t', '--template-pynast', help='Path to a template alignment file if available for Pynast (format: fasta).' )
    group_input.add_argument( '-b', '--biomfile', help='Path to the abundance biom file.' )
        
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-o','--out-tree', default='tree.nwk', help="Path to store resulting Newick tree file. [Default: %(default)s]" )
    group_output.add_argument('-s','--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )    
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)
    
    # Temporary files
    temps=TmpFiles(os.path.split(args.out_tree)[0])
    filename_prefix = ".".join(os.path.split(args.input_otu)[1].split('.')[:-1])
    if args.template_pynast is None:
        align= os.path.join(temps.tmp_dir ,filename_prefix+ '_mafft_aligned.fasta')   
        temps.files.append(align)
        pynast_fail=None
    else:
        align = os.path.join(temps.tmp_dir , filename_prefix+ '_pynast_aligned.fasta')
        pynast_fail = os.path.join(temps.tmp_dir , filename_prefix+'_pynast_fail.fasta')
        pynast_log = os.path.join(temps.tmp_dir , filename_prefix+'_pynast_log.txt') 
        temps.files.append(align)
        temps.files.append(pynast_fail)
        temps.files.append(pynast_log)		
    
    # Process 
    try:        
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        nb_seq = get_fasta_nb_seq(args.input_otu)
        Logger.static_write(args.log_file, "Number of input OTUs sequences: " + str(nb_seq) + "\n\n")
        if nb_seq >10000:
            raise Exception( "FROGS Tree is only working on less than 10 000 sequences!" )
        if args.template_pynast is None:
            mafftMet=get_methods_mafft(args.input_otu)
            Mafft(mafftMet, args.input_otu, align, args.nb_cpus).submit( args.log_file )
        else:
            min_len=compute_min_sequence_length(args.input_otu)
            Pynast(args.input_otu, args.template_pynast, min_len, align, pynast_fail, pynast_log).submit( args.log_file )
        FastTree(align, args.out_tree).submit( args.log_file )
        write_summary( args.html, pynast_fail, args.biomfile, args.out_tree)
    finally:
        if not args.debug:
            temps.deleteAll()
