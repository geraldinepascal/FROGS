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

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("phyloseq_beta_diversity"):
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
    @summary: Launch Rmarkdown script to present the data beta diversity with phyloseq.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    @return: html file containing the plots
             beta divesity distance matrix tsv file(s)
    """
    def __init__(self, html, phyloseq, varExp, methods, outdir, rmd_stderr):
        """
        @param html: [str] path to store resulting html file.
        @param phyloseq: [str] path to phyloseq object in RData file, the result of FROGS Phyloseq Import Data.
        @param varExp: [str] Experiment variable to split plot.
        @param methods: [str] one or more of beta diversity method.        
        @param outdir: [str] The path to store resulting beta diversity distance matrix.
        @param rmd_stderr: [str] Path to temporary Rmarkdown stderr output file
        """ 
        rmd = os.path.join(CURRENT_DIR, "phyloseq_beta_diversity.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render(' + "'" + rmd + "',knit_root_dir='" + outdir + "',output_file='" + html + \
                       "', params=list(phyloseq='" + phyloseq + "', varExp='" + varExp + "', methods='" + methods + "', libdir ='" + LIBR_DIR + "', version='"+ str(__version__) + "'), intermediates_dir='" + os.path.dirname(html) + "')" + '" 2> ' + rmd_stderr,
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
    parser = argparse.ArgumentParser( description='To present the data beta diversity with phyloseq.')    
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]" )   
    
    parser.add_argument('--varExp', type=str, required=True, default=None, help='The experiment variable you want to analyse. [Default: %(default)s]')
    parser.add_argument('--distance-methods', required=True, type=str, default='bray,cc,unifrac,wunifrac', help='Comma separated values beta diversity methods available in Phyloseq (see https://www.bioconductor.org/packages/devel/bioc/manuals/phyloseq/man/phyloseq.pdf). [Default: %(default)s].')
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--rdata', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data. [Default: %(default)s]" )
    # output
    group_output = parser.add_argument_group( 'Outputs' )    
    group_output.add_argument('--matrix-outdir', required=True, action="store", type=str, help="Path to output matrix file")       
    group_output.add_argument('--html', default='phyloseq_beta_diversity.nb.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    # check parameter
    list_distance=["unifrac","wunifrac","bray","cc","dpcoa","jsd","manhattan","euclidean","canberra","kulczynski","jaccard","gower","altGower","morisita","horn","mountford","raup","binomial","chao","cao","wt","-1","c","wb","rt","I","e","t","me","j","sor","m","-2","co","g","-3","l","19","hk","rlb","sim","gl","z","maximum","binary","minkowski","ANY"]
    
    methods = args.distance_methods.strip() if not args.distance_methods.strip()[-1]=="," else args.distance_methods.strip()[:-1]
    for method in methods.split(","):
        if method not in list_distance:
            raise_exception( Exception( '\n\n#ERROR : Your method "'+str(method)+'", name is not correct !!! Please make sure that it is in the list:'+str(list_distance)+"\n\n"))

    # Process 
    outdir = os.path.abspath(args.matrix_outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    phyloseq=os.path.abspath(args.rdata)
    html=os.path.abspath(args.html)
    try:
        tmpFiles = TmpFiles(os.path.dirname(html))
        rmd_stderr = tmpFiles.add("rmarkdown.stderr")
        Rscript(html, phyloseq, args.varExp, methods, outdir, rmd_stderr).submit( args.log_file )
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
