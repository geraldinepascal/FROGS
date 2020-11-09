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
if CURRENT_DIR.endswith("phyloseq_import"):
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
from frogsBiom import *
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
    def __init__(self, biomfile, samplefile, treefile, html, normalization, phyloseq, ranks, rmd_stderr):
        """
        @param biomfile: [str] The biom file contains the  OTU's informations: abundance and taxonomy. These file is the result of FROGS.
        @param samplefile: [str] The tabular file contains the samples's informations. 
                                 Advice: within SampleID or without SampleID
        @param treefile: [str] The Newick file contains the tree's informations from Frogs Tree.
        @param html: [str] The path to store resulting html file.
        @param normalization: [str] To normalize data before analysis.
        @param phyloseq: [str] The path to store one phyloseq-class object in Rdata file.
        @param ranks: [str] The ordered taxonomic ranks levels stored in BIOM. Each rank is separated by one space.
        @param rmd_stderr: [str] Path to temporary Rmarkdown stderr output file
        """ 
        # rmd = os.path.join(CURRENT_DIR, "r_import_data_notebook.Rmd")
        rmd = os.path.join(CURRENT_DIR, "phyloseq_import_data.Rmd")
        Cmd.__init__( self,
                      'Rscript',
                      'Run r_import_data.Rmd',
                      '-e "rmarkdown::render(' + "'" + rmd + "',output_file='" + html + \
                      "', params=list(biomfile='" + biomfile + "', samplefile='" + samplefile + "', treefile='"+ treefile + \
                      "', normalization=" + normalization + ", outputRdata='" + phyloseq + "', ranks='" + ranks +"', libdir ='" + LIBR_DIR + "'), intermediates_dir='" + os.path.dirname(html) + "')" + '" 2> ' + rmd_stderr,
                       "-e '(sessionInfo()[[1]][13])[[1]][1]; paste(\"Rmarkdown version: \",packageVersion(\"rmarkdown\")) ; library(phyloseq); paste(\"Phyloseq version: \",packageVersion(\"phyloseq\"))'")
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')

class FROGSBiomToStdBiom(Cmd):
    """
    @summary : standardize FROGS biom file by keeping blast affiliation consensus as final taxonomy
    """
    def __init__(self, biom_in, biom_out, blast_metadata):
        """
        @param biom_in : [str] The FROGS input biom file path
        @param biom_out : [str] The output standard biom file path
        @param blast_metadata : [str] Output tsv file containing blast detailed affiliations
        """    
        Cmd.__init__(self,
                    'biom_to_stdBiom.py',
                    'standardize FROGS biom file',
                    '--input-biom ' + biom_in + ' --output-biom ' + biom_out + ' --output-metadata ' + blast_metadata,
                    '--version'
        )

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
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )   
    parser.add_argument( '--version', action='version', version=__version__ )
    parser.add_argument( '-n','--normalization', default=False, action='store_true', help='To normalize data before analysis. Use this option if you didnt do it in FROGS Abundance normalisation. [Default: %(default)s]')
    parser.add_argument( '-r','--ranks', type=str, nargs='*', default=['Kingdom', 'Phylum', 'Class', 'Order','Family','Genus', 'Species'], help='The ordered taxonomic ranks levels stored in BIOM. Each rank is separated by one space. [Default: %(default)s]')      
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-b', '--biomfile', required=True, help='path to the abundance BIOM file.' )
    group_input.add_argument( '-s', '--samplefile', required=True, help='path to sample file (format: TSV).' )
    group_input.add_argument( '-t', '--treefile', default=None, help='path to tree file from FROGS Tree (format: Newick "nhx" or "nwk" ).' )
   
    # output
    group_output = parser.add_argument_group( 'Outputs' ) 
    group_output.add_argument('--rdata', default='phyloseq_data.Rdata', help="path to store phyloseq-class object in Rdata file. [Default: %(default)s]" )
    group_output.add_argument('-o','--html', default='phyloseq_import_summary.nb.html', help="The HTML file containing the graphs. [Default: %(default)s]" )
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several informations on executed commands.')   
    args = parser.parse_args()
    prevent_shell_injections(args)
   
    # Process  
    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
    html=os.path.abspath(args.html)
    phyloseq=os.path.abspath(args.rdata)
    biomfile=os.path.abspath(args.biomfile)
    samplefile=os.path.abspath(args.samplefile)

    biom = BiomIO.from_json(biomfile)
    # biom file need to be standardize
    to_standardize = False
    if not biom.has_metadata("taxonomy") :
        if not biom.has_metadata("blast_taxonomy"):
            raise Exception("\n\n#ERROR : Your biom input file is not comming from FROGS and has no standard taxonomy metadata.\n\n")
        else:
            to_standardize=True
    # check sample names compatibility between input biom and sample metadata file
    sample_metadata_list = set()
    FH_in = open(args.samplefile)
    FH_in.readline()
    for line in FH_in:
        sample_metadata_list.add(line.split()[0])

    biom_sample_list = set([name for name in biom.get_samples_names()])
    sample_metadata_spec = sample_metadata_list.difference(biom_sample_list) 
    sample_biom_spec = biom_sample_list.difference(sample_metadata_list)
    if len(sample_biom_spec) > 0 :
        Logger.static_write(args.log_file, "# WARNING : " + str(len(sample_biom_spec)) + " samples from your biom file are not present in your sample metadata file. They will be excluded from further analysis \n\t" + "; ".join(sample_biom_spec) + "\n\n")
    if len(sample_metadata_spec) > 0 :
       raise Exception( "\n\n#ERROR : " + str(len(sample_metadata_spec)) + " among " + str(len(sample_metadata_list)) + " samples from your sample metadata file are not present in your biom file:\n\t" + ";".join(sample_metadata_spec) + "\nPlease give a sample metadata file that fits your abundance biom file\n\n")

    if (args.treefile is None) :
        treefile="None"
    else:
        treefile=os.path.abspath(args.treefile)
    ranks=" ".join(args.ranks)

    try : 
        tmpFiles = TmpFiles(os.path.dirname(html))
        if to_standardize:
            std_biom = tmpFiles.add(os.path.basename(os.path.splitext(biomfile)[0])+".stdBiom")
            blast_metadata = tmpFiles.add(os.path.basename(os.path.splitext(biomfile)[0])+".blast_metadata")
            FROGSBiomToStdBiom(biomfile, std_biom, blast_metadata).submit(args.log_file)
            biomfile = std_biom
        rmd_stderr = tmpFiles.add("rmarkdown.stderr")
        Rscript(biomfile, samplefile, treefile, html, str(args.normalization).upper(), phyloseq, ranks, rmd_stderr).submit(args.log_file)
    finally :
        if not args.debug:
            tmpFiles.deleteAll()
