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
__author__ = 'Ta Thi Ngan & Maria Bernard INRA - SIGENAE'
__copyright__ = 'Copyright (C) 2017 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH']
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

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
    def __init__(self, html, data, varExp, methods, outdir):
        """
        @param html: [str] path to store resulting html file.
        @param data: [str] path to phyloseq object in RData file, the result of FROGS Phyloseq Import Data.
        @param varExp: [str] Experiment variable to split plot.
        @param methods: [str] one or more of beta diversity method.        
        @param outdir: [str] The path to store resulting beta diversity distance matrix.
        """ 
        rmd = os.path.join(CURRENT_DIR, "r_beta_diversity.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render('+"'"+rmd+"',knit_root_dir='"+outdir+ "',output_file='"+html+"', params=list(data='"+data+"', varExp='"+varExp+ "', methods='"+methods+ "'), intermediates_dir='"+os.path.dirname(html)+"')"+'" 2> /dev/null',
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
    parser.add_argument('--output_dir', required=True, action="store", type=str, help="Path to output matrix file")       
    parser.add_argument('-v', '--varExp', type=str, required=True, default=None, help='The experiment variable you want to analyse.')
    parser.add_argument('-m', '--distance-methods', required=True, type=str, default='bray,cc,unifrac,wunifrac', help='The methods of beta diversity in Phyloseq distance list. Each method is separated by one comma. [Default: %(default)s].')
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-d','--data', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data" )
    # output
    group_output = parser.add_argument_group( 'Outputs' )    
    group_output.add_argument('-o','--html', default='beta_diversity.html', help="Path to store resulting html file. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log_file', default=sys.stdout, help='This output file will contain several information on executed commands.')    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    # check parameter
    if args.distance_methods=="None":
        raise Exception("You must choose at least one method name.")     
    list_distance=["unifrac","wunifrac","bray","cc","dpcoa","jsd","manhattan","euclidean","canberra","kulczynski","jaccard","gower","altGower","morisita","horn","mountford","raup","binomial","chao","cao","w","-1","c","wb","r","I","e","t","me","j","sor","m","-2","co","g","-3","l","19","hk","rlb","sim","gl","z","maximum","binary","minkowski","ANY"]
    
    methods=args.distance_methods.split(",") if not args.distance_methods[-1] == "," else args.distance_methods[:-1].split(",")
    for method in methods:
        if method not in list_distance:
            raise Exception( 'Your method "'+str(method)+'", name is not correct !!! Please make sure that it is in the list:'+str(list_distance))

    outdir = args.output_dir
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    outdir=os.path.abspath(outdir)

    # Process 
    data=os.path.abspath(args.data)
    html=os.path.abspath(args.html)
    Rscript(html, data, args.varExp, args.distance_methods, outdir).submit( args.log_file )
