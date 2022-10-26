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
__version__ = '4.0.1'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("phylo_beta_diversity"):
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
        @param analysis: [str] OTU or FUNC: Type of analysis to be done.
        @param data : [str] OTU: The path of one phyloseq-class object in Rdata file.
        @param var: [str] Experimental variable suspected to have an impact on OTUs abundances.
        @param function_table: [str] FUNC: Path to function prediction abundances table from FROGSFUNC function step.
        @param samplefile: [str]: FUNC: Path to metadata samplefile.
        @param out  : [str] Path to Rdata file storing DESeq2 prepreocessing step.
        @param stderr  : [str] Path to stderr output file
        """ 
        rcode = os.path.join(BIN_DIR, "deseq2_preprocess.R")
        if analysis == "OTU":
            opt = ' --inRdata ' + data
        elif analysis == "FUNC":
            opt = ' --inputFunction ' + function_table + ' --samplefile ' + samplefile

        Cmd.__init__( self,
                      'deseq2_preprocess.R',
                      'Construc DESeq2 object from a Phyloseq one.',
                      ' --analysis ' + analysis + ' --var ' + var + ' --outRdata ' + out + opt + ' 2> ' + stderr,
                      '--version')       
                       
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return : [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

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
    group_input.add_argument('-a', '--analysis', required=True, choices=['OTU', 'FUNC'], help='Type of data to perform the differential analysis. OTU: DESeq2 is run on the OTUs abundances table. FUNC: DESeq2 is run on FROGSFUNC function abundances table (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step).')

    group_input_otu_table = parser.add_argument_group( ' OTU ' )
    group_input_otu_table.add_argument('-d','--data', default=None, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data. Required.")
	
    group_input_function_table = parser.add_argument_group( ' FUNC ' )
    group_input_function_table.add_argument('-f', '--input-function', default=None, help='Input file of metagenome function prediction abundances (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step). Required. (default: %(default)s).')
    group_input_function_table.add_argument('-s', '--samplefile', default=None, help='path to sample file (format: TSV). Required.' )
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--out-Rdata', default='DESeq2_preprocess.Rdata', help="The path to store resulting dataframe of DESeq2. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

	# Check for OTU input
    data = args.data
    if args.analysis == "OTU" and data is None:
        parser.error("\n\n#ERROR : --data is required for OTUs analysis. ")
    elif args.analysis == "OTU":
        data=os.path.abspath(args.data)

	# Check for ITS or 18S input
    elif args.analysis == "FUNC" and (args.input_function is None or args.samplefile is None):
        parser.error("\n\n#ERROR : --input-function-table and --samplefile both required for FROGSFUNC analysis.\n\n")
   
    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    out_Rdata=os.path.abspath(args.out_Rdata)

    tmpFiles = TmpFiles(os.path.dirname(out_Rdata))

    try:
        R_stderr = tmpFiles.add("R.stderr")
        Rscript(args.analysis, data, args.var, args.input_function, args.samplefile, out_Rdata,R_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()

