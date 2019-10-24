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
#
__author__ = ' Ta Thi Ngan & Maria Bernard SIGENAE / Mahendra Mariadassou plateforme Migale '
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.1'
__email__ = 'frogs-support@inra.fr'
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
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

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
    def __init__(self, data, var, out, stderr ):
        """
        @param data : [str] The path of one phyloseq-class object in Rdata file.
        @param model: [str] Experimental variable suspected to have an impact on OTUs abundances.
        @param out  : [str] Path to Rdata file storing DESeq2 prepreocessing step.
        @param stderr  : [str] Path to stderr output file
        """ 
        rcode = os.path.join(CURRENT_DIR, "deseq2_preprocess.R")
        Cmd.__init__( self,
                      'Rscript',
                      'Run deseq2_preprocess.R',
                      '--vanilla ' + rcode + ' --inRdata ' + data + ' --var ' + var + ' --outRdata ' + out + ' 2> ' + stderr,
                      "-e '(sessionInfo()[[1]][13])[[1]][1]; library(DESeq2); paste(\"DESeq2 version: \",packageVersion(\"DESeq2\")); library(phyloseq); paste(\"Phyloseq version: \",packageVersion(\"phyloseq\"))'")       
                       
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
    parser.add_argument('-v', '--var', type=str, required=True, help='Experimental variable suspected to have an impact on OTUs abundances.' )   
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-d','--data', required=True, default=None, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data")

    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--out-Rdata', default='DESeq2_preprocess.Rdata', help="The path to store resulting dataframe of DESeq2. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)
   
    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    data=os.path.abspath(args.data)
    out_Rdata=os.path.abspath(args.out_Rdata)

    tmpFiles = TmpFiles(os.path.dirname(out_Rdata))

    try:
        R_stderr = tmpFiles.add("R.stderr")
        Rscript(data, args.var, out_Rdata,R_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()

