#!/usr/bin/env python3

__author__ = ' Ta Thi Ngan - SIGENAE/GABI & Maria Bernard - SIGENAE/GABI & Mahendra Mariadassou - MaIAGE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse
import pandas as pd

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("deseq2_preprocess"):
    FROGS_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
else:
    FROGS_DIR = os.path.dirname(CURRENT_DIR)

# PATH
BIN_DIR = os.path.abspath(os.path.join(FROGS_DIR, "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']

# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(FROGS_DIR, "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

from frogsUtils import *
from frogsBiom import *
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Rscript(Cmd):
    """
    @summary: Launch Rscript to calcul data frame of DESEq2 from a phyloseq object in RData file, the result of FROGS Phyloseq Import Data.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    """
    def __init__(self, analysis, data, var, function_table, samplefile, out, stderr ):
        """
        @param analysis: [str] ASV or FUNCTION: Type of analysis to be done.
        @param data : [str] [ASV]: The path of one phyloseq-class object in Rdata file.
        @param var: [str] Experimental variable suspected to have an impact on ASVs/FUNCTIONs abundances.
        @param function_table: [str] [FUNCTION]: Path to function prediction abundances table from FROGSFUNC function step.
        @param samplefile: [str]: [FUNCTION]: Path to metadata samplefile.
        @param out  : [str] Path to Rdata file storing DESeq2 prepreocessing step.
        @param stderr  : [str] Path to stderr output file
        """ 
        rcode = os.path.join(BIN_DIR, "deseq2_preprocess.R")
        if analysis == "ASV":
            opt = ' --inRdata ' + data
        elif analysis == "FUNCTION":
            opt = ' --inputFunction ' + function_table + ' --samplefile ' + samplefile

        Cmd.__init__( self,
                      'deseq2_preprocess.R',
                      'Construct DESeq2 object from a Phyloseq one.',
                      ' --analysis ' + analysis + ' --var ' + var + ' --outRdata ' + out + opt + ' 2> ' + stderr,
                      '--version')       
                       
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return : [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

class Tsv2biom(Cmd):
    """
    @summary: Create a temporary biom file for FUNCTION phyloseq data object.
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

class PhyloseqImport(Cmd):
    """
    @summary: import data from two files: biomfile and samplefile into a phyloseq object for FUNCTION analysis.
    """
    def __init__(self, biom_file, sample_file, ranks, out_rdata, out_html, log):
        """
        @param biom_file: [str] Path to biom file of function abundances from frogsfunc_functions.py step.
        @param sample_file: [str] Path to samplefile of metadata.
        @param out_rdata: [str] Phyloseq rdata output object.
        @param log: [str] log file.
        """

        Cmd.__init__(self,
                 'phyloseq_import_data.py',
                 'create phyloseq object like with function abundances and annotation', 
                 ' --input-biom ' + biom_file + ' --sample-metadata-tsv ' + sample_file + ' --ranks ' + ranks + ' --output-phyloseq-rdata ' + out_rdata + ' --html ' + out_html + '  2>> ' + log,
                "--version")

    def get_version(self):
        return Cmd.get_version(self, 'stdout').strip()

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def formate_abundances_file( in_tsv, out_tsv):
    df = pd.read_csv(in_tsv, sep='\t')
    df = df.drop('db_link', axis=1)
    df = df.rename(columns={'classification': '#taxonomy'})
    headers = ['#taxonomy', 'observation_name', 'observation_sum']
    for column in df:
        if column not in headers:
            df[column] = df[column].round(0).astype(int)
    df.to_csv(out_tsv, sep="\t", index=False)

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
   
    # Manage parameters
    parser = argparse.ArgumentParser( description='Launch Rscript to generate dataframe of DESEq2 from a phyloseq object in RData file')
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument('--var-exp', type=str, required=True, help='Experimental variable suspected to have an impact on abundances. \
        You may precise complexe string such as variables with confounding effect (ex: Treatment+Gender or Treatmet*Gender)' )   
    parser.add_argument('--analysis-type', required=True, choices=['ASV', 'FUNCTION'], help='Differential analysis on ASV (see phyloseq_import.py) or on Function abundance (see frogsfunc_functions.py).')

    # Inputs
    group_input = parser.add_argument_group()
    group_input_asv = parser.add_argument_group( '# Inputs for ASV analysis type ' )
    group_input_asv.add_argument('--phyloseq-rdata', help="The path of RData file containing a phyloseq object-the result of phyloseq_import.py. Required." )

    group_input_function = parser.add_argument_group( '# Inputs for FUNCTION analysis type ' )
    group_input_function.add_argument('--input-functions-abund', help='Input file of metagenome function prediction abundances (frogsfunc_functions_unstrat.tsv from frogsfunc_functions.py). Required.')
    group_input_function.add_argument('--sample-metadata-tsv', help='path to sample file (format: TSV). Required.' )
    
    # output
    group_output = parser.add_argument_group( '# Outputs' )

    group_output_fun = parser.add_argument_group( '  ## Outputs specific of FUNCTION analysis type ' )
    group_output_fun.add_argument('--output-phyloseq-rdata', default='phyloseq_fun.Rdata', help="Rdata file path to store phyloseq-class object based on functions abundances and annotation. [Default: %(default)s]" )
    
    group_output.add_argument('--output-deseq-rdata', default=None, help="The path to store resulting dataframe of DESeq2. [Default: %(default)s]" )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several information on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)
    
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

    # Check for ASV input
    data = args.phyloseq_rdata
    if args.analysis_type == "ASV" and data is None:
        parser.error("\n\n#ERROR : --data is required for ASVs analysis. ")
    elif args.analysis_type == "ASV":
        data=os.path.abspath(args.phyloseq_rdata)

    # Check for FUNCTION input
    if args.analysis_type == "FUNCTION":
        if args.input_functions_abund is None or args.sample_metadata_tsv is None:
            parser.error("\n\n#ERROR : --input-functions and --samplefile both required for FROGSFUNC analysis.\n\n")

    # Adapt default output file name
    if args.output_deseq_rdata is None:
        if args.analysis_type == "ASV":
            args.output_deseq_rdata = "asv_dds.Rdata"
        elif args.analysis_type == "FUNCTION":
            args.output_deseq_rdata = "function_dds.Rdata"

    out_Rdata=os.path.abspath(args.output_deseq_rdata)
    tmpFiles = TmpFiles(os.path.dirname(out_Rdata))

    # FUNCTION : phyloseq object generation
    if args.analysis_type == "FUNCTION":
        tmp_function_abund_tostd = tmpFiles.add( "functions_unstrat_toStdbiom.tsv")
        formate_abundances_file(args.input_functions_abund, tmp_function_abund_tostd)

        tmp_function_abundances_biom = tmpFiles.add( "function_abundances.biom")
        Tsv2biom(tmp_function_abund_tostd, tmp_function_abundances_biom).submit( args.log_file)

        # check sample names compatibility between input biom and sample metadata file
        #       - if more samples in abundance file than in sample_metadata ==> ok supplementary samples will be excluded
        #       - if more samples in sample_metadata ==> Error 
        #       - this is the default behavior of phyloseq (and check in phyloseq_import)
        
        sample_metadata_list = set()
        FH_in = open(args.sample_metadata_tsv)
        FH_in.readline()
        for line in FH_in:
            sample_metadata_list.add(line.split()[0]) 

        biom = BiomIO.from_json(tmp_function_abundances_biom)
        biom_sample_list = set([name for name in biom.get_samples_names()])
        sample_metadata_spec = sample_metadata_list.difference(biom_sample_list) 
        sample_biom_spec = biom_sample_list.difference(sample_metadata_list)
        if len(sample_biom_spec) > 0 :
            Logger.static_write(args.log_file, "# WARNING : " + str(len(sample_biom_spec)) + " samples from your biom file are not present in your sample metadata file. They will be excluded from further analysis \n\t" + "; ".join(sample_biom_spec) + "\n\n")
        if len(sample_metadata_spec) > 0 :
           raise_exception( Exception( "\n\n#ERROR : " + str(len(sample_metadata_spec)) + " among " + str(len(sample_metadata_list)) + " samples from your sample metadata file are not present in your biom file:\n\t" + "; ".join(sample_metadata_spec) + "\nPlease give a sample metadata file that fits your abundance biom file\n\n"))

        ranks = " ".join(['Level_4', 'Level_3', 'Level_2', 'Level_1'])
        phyloseq_log = tmpFiles.add( "phyloseq_import.log")
        phyloseq_html = tmpFiles.add( "phyloseq_import.nb.html")
        PhyloseqImport(tmp_function_abundances_biom, args.sample_metadata_tsv, ranks, args.output_phyloseq_rdata, phyloseq_html, phyloseq_log).submit( args.log_file)

    try:
        R_stderr = tmpFiles.add("R.stderr")
        Rscript(args.analysis_type, data, args.var_exp, args.input_functions_abund, args.sample_metadata_tsv, out_Rdata, R_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
