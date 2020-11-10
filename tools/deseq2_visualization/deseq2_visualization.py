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
__author__ = ' Ta Thi Ngan SIGENAE / Mahendra Mariadassou plateforme MIGALE'
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs@toulouse.inrae.fr'
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
LIBR_DIR = os.path.join(LIB_DIR,"external-lib")

from frogsUtils import *
##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################

class Rscript(Cmd):
    """
    @summary: Launch Rscript to visualize Differential gene expression analysis.
    @see:  http://rmarkdown.rstudio.com/
           https://joey711.github.io/phyloseq/
    """
    def __init__(self, phyloseq_data, dds, var, mod1, mod2, padj, html, err ):
        """
        @param data: [str] The path of one phyloseq-class object in Rdata file. 
        @param dds : [str] The path of rdata file containing the DESeqDataSet.
        @param var : [str] The experiment variable.
        @param mod1: [str] one variance of variable that you want to test.
        @param mod2: [str] one other variance of variable that you want to test.
        @param padj: [str] the adjusted p-value.
        @param html: [str] path to store resulting html file.
        @param err:  [str] path to store RScript stderr output 

        """ 
        rmd = os.path.join(CURRENT_DIR, "deseq2_visualization.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run deseq2_visualization.Rmd',
                       '-e "rmarkdown::render(' + "'" + rmd + "', output_file='" + html + "', params=list(phyloseq_data='" + phyloseq_data + "', dds='" + dds + "', var='" + var+"', mod1='" + mod1 + "', mod2='" + mod2 + "', padj=" + str(padj) + "))" + '" 2> ' + err ,
                      "-e '(sessionInfo()[[1]][13])[[1]][1]; library(DESeq2); paste(\"DESeq2 version: \",packageVersion(\"DESeq2\"))'")
                      
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
   
    # Manage parameters
    parser = argparse.ArgumentParser( description='Launch Rmarkdown to visualize differential abundance analysis.')
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument('-v', '--varExp', type=str, required=True, help='variable that you want to test.' )
    parser.add_argument('-m1', '--mod1', type=str, default="None", help='one value of the tested variable you want to compare (if more than 2 value in your experiement variable analyzed.)' )
    parser.add_argument('-m2', '--mod2', type=str, default="None", help='second value of the tested variable you want to compare.(if more than 2 value in your experiement variable analyzed.)' )
    parser.add_argument('-pa', '--padj', type=float, default=0.05, help='the adjusted p-value threshold to defined OTU as differentially abundant. [Default: %(default)s]' )

    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-p','--phyloseqData', required=True, help="The path to the RData file containing a phyloseq object (result of FROGS Phyloseq Import Data)")
    group_input.add_argument('-d','--dds', required=True, help="The path to the Rdata file containing the DESeq dds object (result of FROGS DESeq2 Preprocess)")   

    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--html', default='DESeq2_visualization.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)
   
    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    phyloseq_data=os.path.abspath(args.phyloseqData)
    dds=os.path.abspath(args.dds)
    html=os.path.abspath(args.html)
    
    tmpFiles = TmpFiles(os.path.dirname(html))

    try:
        R_stderr = tmpFiles.add("R.stderr")
        Rscript(phyloseq_data, dds, args.varExp, args.mod1, args.mod2, args.padj, html, R_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
