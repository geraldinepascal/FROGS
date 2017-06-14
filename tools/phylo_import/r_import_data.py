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
__author__ = ' Ta Thi Ngan & Maria Bernard INRA - SIGENAE '
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
    @summary: Launch Rmarkdown script to import data from three files: biomfile, samplefile, treefile into a phyloseq object.
    @see: http://rmarkdown.rstudio.com/
          https://joey711.github.io/phyloseq/
    """
    def __init__(self, biomfile, samplefile, treefile, html, normalization, data, ranks):
        """
        @param biomfile: [str] The biom file contains the  OTU's informations: abundance and taxonomy. These file is the result of FROGS.
        @param samplefile: [str] The tabular file contains the samples's informations. 
                                 Advice: within SampleID or without SampleID
        @param treefile: [str] The Newick file contains the tree's informations from Frogs Tree.
        @param html: [str] The path to store resulting html file.
        @param normalization: [str] To normalize data before analysis.
        @param data: [str] The path to store one phyloseq-class object in Rdata file.
        @param ranks: [str] The ordered taxonomic ranks levels stored in BIOM. Each rank is separated by one space.
        """ 
        rmd = os.path.join(CURRENT_DIR, "r_import_data.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run r_import_data.Rmd',
                      '-e "rmarkdown::render('+"'"+rmd+"',output_file='"+html+"', params=list(biomfile='"+biomfile+"', samplefile='"+samplefile+"', treefile='"+treefile+"', normalization="+normalization+", outputRdata='"+data+"', ranks='"+ranks+"'))"+'" 2> /dev/null',
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
    parser = argparse.ArgumentParser( description='Launch Rmardown script to import data from 3 files: biomfile, samplefile, treefile into a phyloseq object')
    parser.add_argument( '-n','--normalization', default=False, action='store_true', help='To normalize data before analysis. Use this option if you didnt do it in FROGS Abundance normalisation. [Default: %(default)s]')
    parser.add_argument( '-r','--ranks', type=str, nargs='*', default=['Kingdom', 'Phylum', 'Class', 'Order','Family','Genus', 'Species'], help='The ordered taxonomic ranks levels stored in BIOM. Each rank is separated by one space. [Default: %(default)s]')      
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-b', '--biomfile', required=True, help='path to biom file (format: biom1). These file is the result of FROGS.' )
    group_input.add_argument( '-s', '--samplefile', required=True, help='path to sample file (format: tabular).' )
    group_input.add_argument( '-t', '--treefile', default=None, help='path to tree file from FROGS Tree (format: Newich "nhx" or "nwk" ).' )
   
    # output
    group_output = parser.add_argument_group( 'Outputs' ) 
    group_output.add_argument('-d','--data', default='phyloseq_data.Rdata', help="path to store phyloseq-class object in Rdata file. [Default: %(default)s]" )
    group_output.add_argument('-o','--html', default='summary.html', help="path to store resulting html file. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log_file', default=sys.stdout, help='This output file will contain several information on executed commands.')   
    args = parser.parse_args()
    prevent_shell_injections(args)
   
    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    html=os.path.abspath(args.html)
    data=os.path.abspath(args.data)
    biomfile=os.path.abspath(args.biomfile)
    samplefile=os.path.abspath(args.samplefile)
    filename_treefile = ".".join(os.path.split(args.treefile)[1].split('.')[:-1])
    if (args.treefile is None) or (args.treefile =='None') or (filename_treefile == ""):
        treefile="None"
    else:
        treefile=os.path.abspath(args.treefile)
    ranks=" ".join(args.ranks)
    Rscript(biomfile, samplefile, treefile, html, str(args.normalization).upper(), data, str(ranks)).submit(args.log_file)
