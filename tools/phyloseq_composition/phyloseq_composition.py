#!/usr/bin/env python3
#
# Copyright (C) 2018 INRA
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
__author__ = ' Ta Thi Ngan & Maria Bernard INRA - SIGENAE '
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("phyloseq_composition"):
    FROGS_DIR = os.path.dirname(os.path.dirname(CURRENT_DIR))
else:
    FROGS_DIR = os.path.dirname(CURRENT_DIR)

# PATH
BIN_DIR = os.path.abspath(os.path.join(FROGS_DIR, "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
APP_DIR = os.path.abspath(os.path.join(FROGS_DIR, "app"))
os.environ['PATH'] = APP_DIR + os.pathsep + os.environ['PATH']
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
    @summary: Launch Rmarkdown script to present the composition of data with package phyloseq
    @see: http://rmarkdown.rstudio.com/ 
          https://joey711.github.io/phyloseq/
    @return: bar plot and composition plot in one html file.
    """
    def __init__(self, html, phyloseq, varExp,taxaRank1,taxaSet1,taxaRank2,numberOfTaxa, rmd_stderr):
        """
        @param html: [str] path to store resulting html file.
        @param phyloseq: [str] path to phyloseq object in RData file, the result of FROGS Phyloseq Import Data.
        @param varExp: [str] The experiment variable used to split plot.
        @param taxaRank1: [str] Taxonomic rank name to subset your data. 
        @param taxaSet1: [str] The taxon name among taxaRank1 to subset your data.
        @param taxaRank2: [str] The sub taxonomic rank name to aggregate your data.
        @param numberOfTaxa: [str] The number of the most abundant taxa to keep at level taxaRank2.
        @param rmd_stderr: [str] Path to temporary Rmarkdown stderr output file
        """ 
        rmd = os.path.join(CURRENT_DIR, "phyloseq_composition.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render(' + "'" + rmd + "',output_file='" + html + \
                       "', params=list(phyloseq='" + phyloseq + "', varExp='" + varExp + \
                       "',taxaRank1='" + taxaRank1 + "',taxaSet1='" + taxaSet1 + "',taxaRank2='" + taxaRank2 + \
                       "',numberOfTaxa=" + str(numberOfTaxa) + ", libdir ='" + LIBR_DIR + "'), intermediates_dir='" + os.path.dirname(html) + "')" + '" 2> ' + rmd_stderr,
                       "-e '(sessionInfo()[[1]][13])[[1]][1]; paste(\"Rmarkdown version: \",packageVersion(\"rmarkdown\")) ; library(phyloseq); paste(\"Phyloseq version: \",packageVersion(\"phyloseq\"))'")
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
    parser = argparse.ArgumentParser( description='Present the composition of data with package phyloseq' )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument('-v', '--varExp', type=str, required=True, help='The experiment variable used to split plot.' )
    parser.add_argument('-r1', '--taxaRank1', type=str, required=True, help='Select taxonomic rank name to subset your data. [ex: Kingdom]' )
    parser.add_argument('-s1', '--taxaSet1', type=str ,nargs='*', required=True, help='Select taxon name among taxaRank1 to subset your data. [ex: Bacteria]' )
    parser.add_argument('-r2', '--taxaRank2', type=str,required=True, help='Select sub taxonomic rank name to aggregate your data. [ex: Phylum]" ' )
    parser.add_argument('-n', '--numberOfTaxa', type=int, required=True, help='The number of the most abundant taxa to keep at taxaRank2. [ex: 9]" ' )
  
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-r','--rdata', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data")

    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--html', default='phyloseq_composition.nb.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    # Process 
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    html=os.path.abspath(args.html)
    phyloseq=os.path.abspath(args.rdata)
    taxaSet1=" ".join(args.taxaSet1)

    try : 
        tmpFiles = TmpFiles(os.path.dirname(html))
        rmd_stderr = tmpFiles.add("rmarkdown.stderr")
        Rscript(html, phyloseq, args.varExp, args.taxaRank1.strip(), str(taxaSet1.strip()), args.taxaRank2.strip(), args.numberOfTaxa, rmd_stderr).submit( args.log_file )

    finally :
        if not args.debug:
            tmpFiles.deleteAll()
