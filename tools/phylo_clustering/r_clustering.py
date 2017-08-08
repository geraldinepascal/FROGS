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
__email__ = 'frogs@inra.fr'
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
    @summary: Launch Rmarkdown script to clustering of samples (from beta diversity distances) using different linkage method.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    @return: the html file containing the clustering plot.
    """
    def __init__(self, html, data, varExp, distance):
        """
        @params html: [str] Path to store resulting html file.
        @param data: [str] One phyloseq object in Rdata file, the result of FROGS Phyloseq Import Data.
        @param varExp: [str] The experiment variable.
        @param distance: [str] Path of data file containing beta diversity distance matrix. These file is the result of FROGS Phyloseq Beta Diversity. 
        """ 
        rmd = os.path.join(CURRENT_DIR, "r_clustering.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run 1 code Rmarkdown',
                       '-e "rmarkdown::render('+"'"+rmd+"',output_file='"+html+"', params=list(data='"+data+"', varExp='"+varExp+"',distance='"+distance+"'), intermediates_dir='"+os.path.dirname(html)+"')"+'" 2> /dev/null',
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
    parser = argparse.ArgumentParser( description='Clustering of samples using different linkage method.' )
    parser.add_argument('-v', '--varExp', type=str, required=True, default=None, help='The experiment variable you want to analyse.')
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument('-r','--rdata', required=True, default=None, help="The path of RData file containing a phyloseq object-the result of FROGS Phyloseq Import Data" )
    group_input.add_argument('-d','--distance-matrix', required=True, default=None, help="The path of data file containing beta diversity distance matrix. These file is the result of FROGS Phyloseq Beta Diversity." ) 

    # output
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument('-o','--html', default='clustering.html', help="The path to store resulting html file. [Default: %(default)s]" )   
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')    
    args = parser.parse_args()
    prevent_shell_injections(args)   
    # Process 
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    html=os.path.abspath(args.html)
    data=os.path.abspath(args.rdata)
    distance=os.path.abspath(args.distance_matrix)
    Rscript(html, data, args.varExp, distance).submit( args.log_file )
