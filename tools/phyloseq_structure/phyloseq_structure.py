#!/usr/bin/env python3

__author__ = 'Ta Thi Ngan - SIGENAE/GABI & Maria Bernard - SIGENAE/GABI'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("phyloseq_structure"):
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
    @summary: Launch Rmarkdown script to visualise data structure with phyloseq.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    @return: the html file containing the plots.
    """
    def __init__(self, html, phyloseq, varExp, method, distance, rmd_stderr):
        """
        @params html: [str] Path to store resulting html file containing ordination plot and heatmap plot.
        @param phyloseq: [str] One phyloseq object in rdata file, the result of FROGS Phyloseq Import Data.
        @param varExp: [str] The experiment variable.
        @param method: [str] The ordination method,the method in list (MDS, NMDS, DCA, CCA, RDA, CAP, DPCoA).
        @param distance: [str] Path of data file containing beta diversity distance matrix. These file is the result of FROGS Phyloseq Beta Diversity.
        @param rmd_stderr: [str] Path to temporary Rmarkdown stderr output file
        """ 
        rmd = os.path.join(CURRENT_DIR, "phyloseq_structure.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render(' + "'" + rmd + "',output_file='" + html + \
                       "', params=list(phyloseq='" + phyloseq + "', varExp='" + varExp + "',method='" + method + "',distance='" + distance + "', libdir ='" + LIBR_DIR + "', version='"+ str(__version__) + "'), intermediates_dir='" + os.path.dirname(html) + "')" + '" 2> ' + rmd_stderr,
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
    parser = argparse.ArgumentParser( description='Visulization of data structure with heatmap plot and ordination plot of Phyloseq.')
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )   
    
    parser.add_argument('--varExp', type=str, required=True, default=None, help='The experiment variable you want to analyse. [Default: %(default)s]')
    parser.add_argument('--ordination-method', type=str, default='MDS', choices=["MDS", "NMDS", "DPCoA", "PCoA"], help="The ordination methods. [Default: %(default)s]")
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--rdata', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data. [Default: %(default)s]" )
    group_input.add_argument('--distance-matrix', required=True, default=None, help="Path of data file containing beta diversity distance matrix. These file is the result of FROGS Phyloseq Beta Diversity. [Default: %(default)s]") 
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--html', default='phyloseq_structure.nb.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    # Process 
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    html=os.path.abspath(args.html)
    phyloseq=os.path.abspath(args.rdata)
    distance=os.path.abspath(args.distance_matrix)
    try:
        tmpFiles = TmpFiles(os.path.dirname(html))
        rmd_stderr = tmpFiles.add("rmarkdown.stderr")
        Rscript(html, phyloseq, args.varExp, args.ordination_method, distance,rmd_stderr).submit( args.log_file )
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
