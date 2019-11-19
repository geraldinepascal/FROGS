#!/usr/bin/env python2.7
#
# Copyright (C) 2019 INRA
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

__author__ = 'Maria Bernard - SIGENAE Jouy en Josas'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inra.fr'
__status__ = 'beta'

import sys
import os
import re
import time
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH
os.environ['PATH'] = CURRENT_DIR + os.pathsep + os.environ['PATH']

# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(CURRENT_DIR)), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsBiom import BiomIO
from frogsSequenceIO import FastaIO


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Creates the count file from grinder distribution (the expected BIOM) and compares it with the input BIOM.' )
    parser.add_argument( '-t', '--tmp-dir', default=os.getcwd(), help='The tmp directory.' )
    parser.add_argument( '-s', '--sample-sep', default=os.getcwd(), help='Sample name separator in grinder filename (ex: "_" in "spl1_ranks.tsv").' )
    parser.add_argument( '-k', '--taxonomy-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in checked biom. [Default: taxonomy]" )
    parser.add_argument( '-m', '--multi-affiliations', action='store_true', help='The taxonomy is produced by FROGS multi-affiliations.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '--databank', required=True, help='Path to the databank source for simulated sequence (format: fasta). The description of sequences must be the taxonomy.' ) 
    group_input.add_argument( '--checked-biom', required=True, help='Path to the abundance file produced by assessed workflow on simulated data (format: BIOM).' )
    group_input.add_argument( '--checked-fasta', required=False, help='Path to the sequence file produced by assessed workflow on simulated data (format: fasta).' )
    group_input.add_argument( '--grinder-dir', required=True, help='Path to the directory with count by initial sequence in simulation.' )
    args = parser.parse_args()

    real_biom = os.path.join( args.tmp_dir, str(time.time()) + "_" + str(os.getpid()) + "_real.biom" )
    checked_biom = os.path.join( args.tmp_dir, str(time.time()) + "_" + str(os.getpid()) + "_checked.biom" )
    
    # List samples
    samples = list()
    for dirname, dirnames, filenames in os.walk( args.grinder_dir ):
        for filename in filenames:
            if filename.endswith( "_ranks.tsv" ):
                sample_name = filename.split("_ranks.tsv")[0].split(args.sample_sep)[0]
                samples.append( {'name':sample_name, 'path':os.path.join(dirname, filename)} )
    
    # Grinder to BIOM
    cmd_grinder2biom = "grinder2biom.py" + \
        " --affiliation " + os.path.abspath(args.databank) + \
        " --output " + real_biom + \
        " --samples"
    for current_sample in samples:
        cmd_grinder2biom += " '" + current_sample['name'] + ":" + current_sample['path'] + "'"
    subprocess.check_call( cmd_grinder2biom, shell=True )

    # Global Comparison of tax
    cmd_Globalcompare = "its_GlobalTaxCmp.py" \
            + " --real-biom " + os.path.abspath(real_biom) \
            + " --real-tax-key 'real_taxonomy'" \
            + " --checked-biom " + os.path.abspath(args.checked_biom) \
            + " --checked-tax-key '" + args.taxonomy_key + "'" \
            + (" --multi-affiliations" if args.multi_affiliations else "")
    print subprocess.check_output( cmd_Globalcompare, shell=True )
    print ""
    

    # Compare detailed expected to obtained
    for current_sample in samples:
        print "#Sample "+current_sample['name']
        cmd_compareSample = "its_biomCmpTax.py" \
            + " --real-biom " + os.path.abspath(real_biom) \
            + " --real-tax-key 'real_taxonomy'" \
            + " --checked-biom " + os.path.abspath(args.checked_biom) \
            + " --checked-tax-key '" + args.taxonomy_key + "'" \
            + (" --multi-affiliations" if args.multi_affiliations else "") \
            + " --sample " + current_sample['name']
        #~ print cmd_compareSample
        print subprocess.check_output( cmd_compareSample, shell=True )
        print ""

    # Remove tmp files
    os.remove( real_biom )
