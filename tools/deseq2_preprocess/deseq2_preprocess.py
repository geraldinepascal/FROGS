#!/usr/bin/env python3
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
#
__author__ = ' Ta Thi Ngan & Maria Bernard SIGENAE / Mahendra Mariadassou plateforme Migale '
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '4.1.0'
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

# LIBR
# LIBR_DIR = os.path.join(LIB_DIR,"external-lib")

from frogsUtils import *
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
                 'predict gene copy number per sequence.', 
                 ' -b ' + biom_file + ' -s ' + sample_file + ' --ranks ' + ranks + ' --rdata ' + out_rdata + ' --html ' + out_html + '  2>> ' + log,
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
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument('-v', '--var', type=str, required=True, help='Experimental variable suspected to have an impact on abundances. \
        You may precise complexe string such as variables with confounding effect (ex: Treatment+Gender or Treatmet*Gender)' )   
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-a', '--analysis', required=True, choices=['ASV', 'FUNCTION'], help='Type of data to perform the differential analysis. ASV: DESeq2 is run on the ASVs abundances table. FUNCTION: DESeq2 is run on FROGSFUNC function abundances table (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step).')

    group_input_asv_table = parser.add_argument_group( ' ASV ' )
    group_input_asv_table.add_argument('-d','--data', default=None, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data. Required.")

    group_input_function_table = parser.add_argument_group( ' FUNCTION ' )
    group_input_function_table.add_argument('-f', '--input-functions', default=None, help='Input file of metagenome function prediction abundances (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step). Required. (default: %(default)s).')
    group_input_function_table.add_argument('-s', '--samplefile', default=None, help='path to sample file (format: TSV). Required.' )
    group_input_function_table.add_argument('--out-Phyloseq', default='function_data.Rdata', help="path to store phyloseq-class object in Rdata file. [Default: %(default)s]" )
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--out-Rdata', default=None, help="The path to store resulting dataframe of DESeq2. [Default: %(default)s]" )
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)
    
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

    # Check for ASV input
    data = args.data
    if args.analysis == "ASV" and data is None:
        parser.error("\n\n#ERROR : --data is required for ASVs analysis. ")
    elif args.analysis == "ASV":
        data=os.path.abspath(args.data)

    if args.out_Rdata is None:
        if args.analysis == "ASV":
            args.out_Rdata = "asv_dds.Rdata"
        elif args.analysis == "FUNCTION":
            args.out_Rdata = "function_dds.Rdata"

    out_Rdata=os.path.abspath(args.out_Rdata)
    tmpFiles = TmpFiles(os.path.dirname(out_Rdata))

    # Check for ITS or 18S input
    if args.analysis == "FUNCTION":
        if args.input_functions is None or args.samplefile is None:
            parser.error("\n\n#ERROR : --input-functions and --samplefile both required for FROGSFUNC analysis.\n\n")

        tmp_function_abund_tostd = tmpFiles.add( "functions_unstrat_toStdbiom.tsv")
        formate_abundances_file(args.input_functions, tmp_function_abund_tostd)

        tmp_function_abundances_biom = tmpFiles.add( "function_abundances.biom")
        Tsv2biom(tmp_function_abund_tostd, tmp_function_abundances_biom).submit( args.log_file)

        ranks = " ".join(['Level_4', 'Level_3', 'Level_2', 'Level_1'])
        phyloseq_log = tmpFiles.add( "phyloseq_import.log")
        phyloseq_html = tmpFiles.add( "phyloseq_import.nb.html")
        PhyloseqImport(tmp_function_abundances_biom, args.samplefile, ranks, args.out_Phyloseq, phyloseq_html, phyloseq_log).submit( args.log_file)

    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

    try:
        R_stderr = tmpFiles.add("R.stderr")
        Rscript(args.analysis, data, args.var, args.input_functions, args.samplefile, out_Rdata, R_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
