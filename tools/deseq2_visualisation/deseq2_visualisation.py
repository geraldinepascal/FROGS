#!/usr/bin/env python3

__author__ = 'Ta Thi Ngan - SIGENAE/GABI & Mahendra Mariadassou - MaIAGE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '5.0.0'
__email__ = 'frogs@toulouse.inrae.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
FROGS_DIR=""
if CURRENT_DIR.endswith("deseq2_visualisation"):
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
    @summary: Launch Rscript to visualise Differential gene expression analysis.
    @see:  http://rmarkdown.rstudio.com/
           https://joey711.github.io/phyloseq/
    """
    def __init__(self, abundance_data, dds, var, mod1, mod2, padj, html, analysis, err, ipath_over, ipath_under, svg_file_over, svg_file_under ):
        """
        @param data: [str] The path of ASV or FUNCTION abundance table in Rdata file. 
        @param dds : [str] The path of rdata file containing the DESeqDataSet.
        @param var : [str] The experiment variable.
        @param mod1: [str] one variance of variable that you want to test.
        @param mod2: [str] one other variance of variable that you want to test.
        @param padj: [str] the adjusted p-value.
        @param html: [str] path to store resulting html file.
        @param analysis: [str] Either ASV or FUNCTION.
        @param over_svg: [str] Path to temporary svg file (FUNCTION analysis)
        @param under_svg: [str] Path to temporary svg file (FUNCTION analysis)
        @param err:  [str] path to store RScript stderr output

        """ 
        rmd = os.path.join(CURRENT_DIR, "deseq2_visualisation.Rmd")
        if  analysis == "FUNCTION":
            opt = ", ipath_over='" + ipath_over + "', ipath_under='" + ipath_under + "', svg_file_over='" + svg_file_over + "', svg_file_under='" + svg_file_under + "'"
        else:
            opt = ""
            
        Cmd.__init__( self,
                      'Rscript',
                      'Run deseq2_visualisation.Rmd',
                       '-e "rmarkdown::render(' + "'" + rmd + "', output_file='" + html + "', params=list(abundance_data='" + abundance_data + "', analysis='" + analysis + "', dds='" + dds + "', var='" + var+"', mod1='" + mod1 + "', mod2='" + mod2 + "', padj_th=" + str(padj) + ", libdir ='" + LIB_DIR + "'" + ", version='"+ str(__version__) + "'" + opt + "), intermediates_dir='" + os.path.dirname(html) +"')" + '" 2> ' + err ,
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
    parser = argparse.ArgumentParser( description='Launch Rmarkdown to visualise differential abundance analysis.')
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    
    parser.add_argument('--var', type=str, required=True, help='variable that you want to test.' )
    parser.add_argument('--mod1', type=str, default="None", help='one value of the tested variable you want to compare (if more than 2 value in your experiement variable analyzed.) [Default: %(default)s]' )
    parser.add_argument('--mod2', type=str, default="None", help='second value of the tested variable you want to compare.(if more than 2 value in your experiement variable analyzed.) [Default: %(default)s]' )
    parser.add_argument('--padj', type=float, default=0.05, help='the adjusted p-value threshold to defined ASV as differentially abundant. [Default: %(default)s]' )
    parser.add_argument('--analysis', default="ASV", required=True, choices=['ASV', 'FUNCTION'], help='Type of data to perform the differential analysis. ASV: DESeq2 is run on the ASVs abundances table. FUNC: DESeq2 is run on FROGSFUNC function abundances table (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step). [Default: %(default)s]')

    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('--abundanceData', required=True, help="The path to the RData file containing the ASV/FUNCTION abundances table. (result of FROGS Phyloseq Import Data)")
    group_input.add_argument('--dds', required=True, help="The path to the Rdata file containing the DESeq dds object (result of FROGS DESeq2 Preprocess)")   
    
    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('--ipath-over', default=None, help="The tsv file of over abundants functions (FUNCTION analysis only) [Default: %(default)s]" )
    group_output.add_argument('--ipath-under', default=None, help="The tsv file of under abundants functions (FUNCTION analysis only) [Default: %(default)s]" )
    group_output.add_argument('--html', default='DESeq2_visualisation.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument('--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)
    output_dir = os.path.dirname(os.path.abspath(args.html))

    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    abundance_data=os.path.abspath(args.abundanceData)
    dds=os.path.abspath(args.dds)
    html=os.path.abspath(args.html)
    tmpFiles = TmpFiles(os.path.dirname(html))

    try:
        R_stderr = tmpFiles.add("R.stderr")
        html=os.path.abspath(args.html)
        if args.analysis == "ASV":
            if args.ipath_over is not None or args.ipath_under is not None:
                parser.error("\n\n#ERROR : --ipath-over and --ipath-under only available for FUNCTION analysis. ")
            Rscript(abundance_data, dds, args.var, args.mod1, args.mod2, args.padj, html, args.analysis, R_stderr, None, None, None, None ).submit(args.log_file)
        elif args.analysis == "FUNCTION":
            svg_ipath_file_over  = os.path.abspath(output_dir + "/" +  "ipath_over.svg")
            svg_ipath_file_under  = os.path.abspath(output_dir + "/" +  "ipath_under.svg")
            if args.ipath_over is None:
                args.ipath_over = "ipath_over.tsv"
            args.ipath_over =  os.path.abspath(output_dir + "/" + os.path.basename(args.ipath_over))
            if args.ipath_under is None:
                args.ipath_under = "ipath_under.tsv"
            args.ipath_under = os.path.abspath(output_dir + "/" +  os.path.basename(args.ipath_under))

            Rscript(abundance_data, dds, args.var, args.mod1, args.mod2, args.padj, html, args.analysis, R_stderr, args.ipath_over, args.ipath_under, svg_ipath_file_over, svg_ipath_file_under).submit(args.log_file)
    
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
