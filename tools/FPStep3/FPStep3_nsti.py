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
# Pour les parse des dataframe
from numpy import median
import pandas as pd
# Pour deziper
import gzip
# Pour les temp_directory
from tempfile import gettempdir
from shutil import rmtree

#from get_version import get_version

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

#import picrust2
#from picrust2.metagenome_pipeline import run_metagenome_pipeline
#from picrust2.util import check_files_exist
##################################################################################################################################################
#
# COMMAND LINES #  metagenome_pipeline.py -i study_seqs.biom \  -m marker_nsti_predicted.tsv.gz \ -f KO_predicted.tsv.gz \   -o KO_metagenome_out
                
##################################################################################################################################################
class metagenome_pipeline(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: predict gene abudance
    @see: https://github.com/picrust/picrust2/wiki
    """
    
    def __init__(self, metagenomeMet, input_biom, marker, function, max_nsti, out_dir, stdout):
    
        Cmd.__init__(self,
                 'metagenome_pipeline.py ',
                 'predict abundance gene', 
                  "" + metagenomeMet  +" -i "+ str(input_biom)+" -m "+ str(marker) +" -f "+ str(function)+" --max_nsti "+ str(max_nsti)+" -o "+ str(out_dir) +' 2> ' + stdout,
                "--version") 
      
    # get version pose un probléme de tableau
    # def get_version(self):
    #     """
    #     @summary: Returns the program version number.
    #     @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
    #     """

    #     return Cmd.get_version(self, 'stdout').split()[1].strip() 
"""
_Début: 
__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse and Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs-support@inra.fr'
__status__ = 'prod'

"""
class Biom2tsv(Cmd):
    """
    @summary: Converts BIOM file to TSV file.
    @note: taxonomyRDP seedID seedSequence blastSubject blastEvalue blastLength blastPercentCoverage blastPercentIdentity blastTaxonomy OTUname SommeCount sample_count
    """
    def __init__( self, out_tsv, in_biom, input_fasta=None ):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv: [str] Path to output TSV file.
        """
        # Sequence file option
        sequence_file_opt = "" if input_fasta is None else " --input-fasta " + input_fasta

        # Check the metadata
        biom = BiomIO.from_json( in_biom )
        obs = biom.rows[0]
        conversion_tags = ""
        if biom.has_observation_metadata( 'comment' ) :
            conversion_tags += "'comment' "
        if biom.has_observation_metadata( 'rdp_taxonomy' ) and biom.has_observation_metadata( 'rdp_bootstrap' ):
            conversion_tags += "'@rdp_tax_and_bootstrap' "
        if biom.has_observation_metadata( 'blast_taxonomy' ):
            conversion_tags += "'blast_taxonomy' "
        if biom.has_observation_metadata( 'blast_affiliations' ):
            conversion_tags += "'@blast_subject' "
            conversion_tags += "'@blast_perc_identity' "
            conversion_tags += "'@blast_perc_query_coverage' "
            conversion_tags += "'@blast_evalue' "
            conversion_tags += "'@blast_aln_length' "
        if biom.has_observation_metadata( 'seed_id' ):
            conversion_tags += "'seed_id' "
        if input_fasta is not None:
            conversion_tags += "'@seed_sequence' "

        frogs_metadata = ["comment", "rdp_taxonomy", "rdp_bootstrap","blast_taxonomy","blast_affiliations","seed_id"]
        if biom.get_observation_metadata(obs["id"]) != None:
            for metadata in biom.get_observation_metadata(obs["id"]):
                if metadata not in frogs_metadata : 
                    conversion_tags += "'"+metadata+"' "
                
        conversion_tags += "'@observation_name' '@observation_sum' '@sample_count'"

        # Set command
        Cmd.__init__( self,
                      'biom2tsv.py',
                      'Converts a BIOM file in TSV file.',
                      "--input-file " + in_biom + sequence_file_opt + " --output-file " + out_tsv + " --fields " + conversion_tags,
                      '--version' )


class Biom2multiAffi(Cmd):
    """
    @summary: Extracts multi-affiliations from a FROGS BIOM file.
    """
    def __init__( self, out_tsv, in_biom ):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv: [str] Path to output TSV file.
        """
        # Set command
        Cmd.__init__( self,
                      'multiAffiFromBiom.py',
                      'Extracts multi-affiliations data from a FROGS BIOM file.',
                      '--input-file ' + in_biom + ' --output-file ' + out_tsv,
                      '--version' )

"""
_FIN: 
__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse and Maria Bernard - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs-support@inra.fr'
__status__ = 'prod'
"""
##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
"""
_Debut: 
__author__ PICRUST2
"""
def run_metagenome_pipeline(input_seqabun,
                            function,
                            max_nsti,
                            marker=None,
                            min_reads=1,
                            min_samples=1,
                            strat_out=False,
                            wide_table=False,
                            skip_norm=False,
                            out_dir='metagenome_out'):
    '''Main function to run full metagenome pipeline. Meant to run modular
    functions largely listed below. Will return predicted metagenomes
    straitifed and unstratified by contributing genomes (i.e. taxa).'''

    if not marker and not skip_norm:
        sys.exit("Table of predicted marker gene copy numbers is required "
                 "unless --skip_norm is specified.")
    elif marker and skip_norm:
        sys.exit("Table of predicted marker gene copy numbers should not be "
                 "specified when --skip_norm option is set.")

    make_output_dir(out_dir)

    # Initialize empty pandas dataframe to contain NSTI values.
    nsti_val = pd.DataFrame()

    study_seq_counts = read_seqabun(input_seqabun)

    pred_function = pd.read_csv(function, sep="\t", dtype={'sequence': str})
    pred_function.set_index('sequence', drop=True, inplace=True)

    # If NSTI column present then remove all rows with value above specified
    # max value. Also, remove NSTI column (in both dataframes).
    if 'metadata_NSTI' in pred_function.columns:
        pred_function, nsti_val = drop_tips_by_nsti(tab=pred_function,
                                                    nsti_col='metadata_NSTI',
                                                    max_nsti=max_nsti)
    if not skip_norm:
        check_files_exist([marker])
        pred_marker = pd.read_csv(marker, sep="\t", dtype={'sequence': str})
        pred_marker.set_index('sequence', drop=True, inplace=True)

        if 'metadata_NSTI' in pred_marker.columns:
            pred_marker, nsti_val = drop_tips_by_nsti(tab=pred_marker,
                                                      nsti_col='metadata_NSTI',
                                                      max_nsti=max_nsti)

        # Re-order predicted abundance tables to be in same order as study seqs.
        # Also, drop any sequence ids that don't overlap across all dataframes.
        study_seq_counts, pred_function, pred_marker = three_df_index_overlap_sort(study_seq_counts,
                                                                                   pred_function,
                                                                                   pred_marker)
        norm_output = path.join(out_dir, "seqtab_norm.tsv.gz")

        # Normalize input study sequence abundances by predicted abundance of
        # marker genes and output normalized table if specified.
        study_seq_counts = norm_by_marker_copies(input_seq_counts=study_seq_counts,
                                                 input_marker_num=pred_marker,
                                                 norm_filename=norm_output)
    else:
        # Get intersecting rows between input files and sort.
        label_overlap = pred_function.index.intersection(study_seq_counts.index).sort_values()

        if len(label_overlap) == 0:
            sys.exit("No sequence ids overlap between both input files.")

        pred_function = pred_function.reindex(label_overlap)
        study_seq_counts = study_seq_counts.reindex(label_overlap)

    # If NSTI column input then output weighted NSTI values.
    if not nsti_val.empty:
        weighted_nsti_out = path.join(out_dir, "weighted_nsti.tsv.gz")
        calc_weighted_nsti(seq_counts=study_seq_counts,
                           nsti_input=nsti_val,
                           outfile=weighted_nsti_out)

    # Determine which sequences should be in the "RARE" category if stratified
    # table is specified.
    if strat_out:
        rare_seqs = []

        if min_reads != 1 or min_samples != 1:
            rare_seqs = id_rare_seqs(in_counts=study_seq_counts,
                                     min_reads=min_reads,
                                     min_samples=min_samples)

    # Generate and return final tables.
    if not strat_out:
        return(None, unstrat_funcs_only_by_samples(pred_function,
                                                   study_seq_counts))

    elif strat_out and not wide_table:
        return(metagenome_contributions(pred_function, study_seq_counts,
                                        rare_seqs),
               unstrat_funcs_only_by_samples(pred_function, study_seq_counts))

    elif strat_out and wide_table:
        return(strat_funcs_by_samples(pred_function, study_seq_counts,
                                      rare_seqs))

def check_files_exist(filepaths):
    '''Takes in a list of filepaths and checks whether they exist. Will
    throw error describing which files do not exist if applicable.'''

    num_nonexist = 0

    missing_files = []

    for filepath in filepaths:

        if not exists(filepath):
            missing_files += [filepath]
            num_nonexist += 1

    if num_nonexist == 0:
        pass
    elif num_nonexist == 1:
        raise ValueError("This input file was not found: " + missing_files[0])
    elif num_nonexist > 1:
        raise ValueError("These input files were not found: " +
                         ", ".join(missing_files))
"""
_FIN: 
__author__ PICRUST2
"""
# Fonction pour formater le biom
"""
# Recupére la sortie de Biom to tsv
# Supprime les deux premiers colonnes (comment et seed_id) et la colone obs_sum
# Utilise pandas
"""
def tsvParse(out_tsv):

    tsv_parse = "out_tsv.tsv"
    # je crée un dataframe
    df = pd.read_csv(out_tsv, header=0, delimiter='\t')
    # Je suprime les mauvaises colonnes
    al = df.drop(df.columns[[0, 1, 3]], axis=1)
    #print(df.drop(df.columns[[0, 1]], axis=1))
    print(al)
    # J écris dans un nouveau fichier
    al.to_csv(tsv_parse, "\t", index=None)

##### Fonction dezip
def dezip(file1, out_file1):

    #file_zip = gzip.GzipFile("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv.gz", 'rb')
    file_zip = gzip.GzipFile(file1, 'rb')
    s = file_zip.read()
    file_zip.close()

    #file_dezip = open ("/Users/moussa/FROGS_moussa/tools/FPStep3/DonneesPourTest/pred_met1.tsv", 'wb')
    file_dezip = open (out_file1, 'wb')
    file_dezip.write(s)
    file_dezip.close()

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":

    # Manage parameters
    parser = argparse.ArgumentParser( description='Phylogenetic tree reconstruction' )

    parser.add_argument('--strat_out', default=False, action='store_true', help='Output table stratified by sequences as well. By ''default this will be in \"contributional\" format ''(i.e. long-format) unless the \"--wide_table\" ''option is set. The startified outfile is named ''\"metagenome_contrib.tsv.gz\" when in long-format.')

    parser.add_argument('--wide_table', default=False, action='store_true', help='Output wide-format stratified table of metagenome ''predictions when \"--strat_out\" is set. This is ''the deprecated method of generating stratified ''tables since it is extremely memory intensive. The ''startified outfile is named ''\"pred_metagenome_strat.tsv.gz\" when this option ' 'is set.')

    parser.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )

    group_input.add_argument('-i', '--input_biom', required=True, help='Input table of sequence abundances (BIOM, TSV, or ''mothur shared file format).')
    
    group_input.add_argument('-f', '--function', metavar='PATH', type=str, help='Table of predicted gene family copy numbers ''(output of hsp.py).')

    group_input.add_argument('-m', '--marker', metavar='PATH', type=str, help='Table of predicted marker gene copy numbers ''(output of hsp.py - typically for 16S).')

    group_input.add_argument('--max_nsti', metavar='FLOAT', type=float, default=False, required=True, help='Sequences with NSTI values above this value will ' 'be excluded (default: %(default)d).')
   # group_input.add_argument('--max_nsti', metavar='FLOAT', type=float, default=2.0, help='Sequences with NSTI values above this value will ' 'be excluded (default: %(default)d).')

    #@group_input.add_argument( '-t', '--out_tsv', default='abundance.tsv', help='This output file will contain the abundance and metadata (format: TSV). [Default: %(default)s]' )

    #group_input.add_argument( '-b', '--input-biom', required=True, help="The abundance file (format: BIOM)." )

    #group_input.add_argument('--min_reads', metavar='INT', type=int, default=1, help='Minimum number of reads across all samples for ''each input ASV. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')

    group_input.add_argument('--min_samples', metavar='INT', type=int, default=1, help='Minimum number of samples that an ASV needs to be ''identfied within. ASVs below this cut-off will be ''counted as part of the \"RARE\" category in the ''stratified output (default: %(default)d).')

    group_input.add_argument('-e', '--en', metavar='PATH', type=str, help='Table of predicted gene family copy numbers and markers ''(output of hsp.py).')

    group_input.add_argument( '--input-fasta', help='The sequences file (format: fasta). If you use this option the sequences will be add in TSV.' )

    #Outputs
    group_output = parser.add_argument_group( 'Outputs')

    #group_output.add_argument('-o', '--out_dir', metavar='PATH', type=str, default='metagenome_out', help='Output directory for metagenome predictions. ''(default: %(default)s).')

    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

    group_output.add_argument( '-t', '--output-tsv', default='abundance.tsv', help='This output file will contain the abundance and metadata (format: TSV). [Default: %(default)s]' )

    group_output.add_argument('--function_abund', metavar='PATH', default='pred_metagenome_unstrat.tsv.gz', help='Output file for metagenome predictions abundance. ''(default: %(default)s).')

    group_output.add_argument('--seqtab', default='seqtab_norm.tsv.gz', help='This output file will contain abundance normalized. ''(default: %(default)s).')

    group_output.add_argument('--weighted', default='weighted_nsti.tsv.gz', help='This output file will contain the nsti value per sample (format: TSV). [Default: %(default)s]' )

    args = parser.parse_args()

    prevent_shell_injections(args)

    stderr = "FPStep3.stderr"
    metagenomeMet = ""
    ### Il faut changer ici : parser le fichier -f (EC, KO, PF,...) produire un autre fichier qui sera pris ###
    # Process 
    try:     
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        # Lancer Biom2tsv : transforme le biom en tsv
        Biom2tsv( args.output_tsv, args.input_biom, args.input_fasta ).submit( args.log_file )
        print("Biom2tsv FINI")
        # Apelle à la fonction de parse du tsv
        tsvParse(args.output_tsv)
        print("tsvParse FINI")

        ########### Création de chemin #############
        #out_dir = 
        out_dir = os.path.abspath(os.path.dirname(args.function_abund)) + "/" +str(time.time()) + "_" + str(os.getpid())
        # Lancer avec le fichier reaction et le fichier marker
        metagenome_pipeline(metagenomeMet, "out_tsv.tsv", args.marker, args.function, args.max_nsti ,out_dir, stderr).submit( args.log_file )
        print("metagenome_Normal FINI")
        ###### 
        os.system("mv " + out_dir +"/"+ "pred_metagenome_unstrat.tsv.gz " + args.function_abund)

        dezip(args.function_abund, "pred_met1.tsv")

        os.system("mv " + out_dir +"/" + "seqtab_norm.tsv.gz " + args.seqtab)
        dezip(args.seqtab, "seqtab1.tsv")
   
        os.system("mv " + out_dir +"/" + "weighted_nsti.tsv.gz " + args.weighted)
        dezip(args.weighted, "weight1.tsv")

        os.system("rm -rf " + out_dir)
        #  Lancer avec le fichier final 
        #metagenome_pipeline(metagenomeMet, args.input_biom, tab_name[ind_no], tab_name[ind_tree], args.out_dir, stderr).submit( args.log_file )
        #print("metagenome_compliqué FINI")
        #print("okok")
    finally:
        print("Partie Finale ")


