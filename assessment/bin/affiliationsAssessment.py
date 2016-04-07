#!/usr/bin/env python2.7
#
# Copyright (C) 2016 INRA
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

__author__ = 'Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'beta'

import os
import re
import time
import argparse
import subprocess
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
    parser.add_argument( '-s', '--sample-sep', default=os.getcwd(), help='Sample name separator in grinder filename (ex: "_" in "spl1_ranks.txt").' )
    parser.add_argument( '-k', '--taxonomy-key', type=str, default="taxonomy", help="The metadata tag used for store taxonomy in checked biom. [Default: taxonomy]" )
    parser.add_argument( '-m', '--multi-affiliations', action='store_true', help='The taxonomy is produced by FROGS multi-affiliations.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '--databank', required=True, help='Path to the databank source for simulated sequence (format: fasta). The description of sequences must be the taxonomy.' ) 
    group_input.add_argument( '--checked-biom', required=True, help='Path to the abundance file produced by assessed workflow on simulated data (format: BIOM).' )
    group_input.add_argument( '--checked-fasta', required=True, help='Path to the sequence file produced by assessed workflow on simulated data (format: fasta).' )
    group_input.add_argument( '--grinder-dir', required=True, help='Path to the directory with count by initial sequence in simulation.' )
    group_input.add_argument( '-u', '--uniq-groups', default=None, help='Path to the file with by line the list of IDs of initial sequences with the same sequence (format: TSV).' )
    args = parser.parse_args()

    real_biom = os.path.join( args.tmp_dir, str(time.time()) + "_" + str(os.getpid()) + "_real.biom" )
    checked_biom = os.path.join( args.tmp_dir, str(time.time()) + "_" + str(os.getpid()) + "_checked.biom" )

    # List samples
    samples = list()
    for dirname, dirnames, filenames in os.walk( args.grinder_dir ):
        for filename in filenames:
            if filename.endswith( "_ranks.txt" ):
                sample_name = filename.split(args.sample_sep)[0]
                samples.append( {'name':sample_name, 'path':os.path.join(dirname, filename)} )

    # Grinder to BIOM
    cmd_grinder2biom = os.path.join(os.path.dirname(os.path.abspath(__file__)), "grinder2biom.py") + \
        " --affiliation " + os.path.abspath(args.databank) + \
        " --output " + real_biom + \
        " --samples"
    for current_sample in samples:
        cmd_grinder2biom += " '" + current_sample['name'] + ":" + current_sample['path'] + "'"
    subprocess.check_call( cmd_grinder2biom, shell=True )

    # Add reference id in checked BIOM
    biom = BiomIO.from_json( args.checked_biom )
    fasta = FastaIO( args.checked_fasta )
    for record in fasta:
        reference = re.search("reference=([^\s]+)", record.description).group(1)
        biom.add_metadata( record.id, "grinder_source", reference, "observation" )
    fasta.close()
    BiomIO.write( checked_biom, biom )
    del(biom)

    # Compare expected to obtained
    for current_sample in samples:
        print current_sample['name']
        cmd_compareSample = os.path.join(os.path.dirname(os.path.abspath(__file__)), "biomCmpTax.py") \
            + " --real-biom " + os.path.abspath(real_biom) \
            + " --real-tax-key 'real_taxonomy'" \
            + " --checked-biom " + os.path.abspath(checked_biom) \
            + " --checked-tax-key '" + args.taxonomy_key + "'" \
            + (" --multi-affiliations" if args.multi_affiliations else "") \
            + (" --uniq-groups " + args.uniq_groups if args.uniq_groups is not None else "") \
            + " --sample " + current_sample['name']
        print subprocess.check_output( cmd_compareSample, shell=True )
        print ""

    # Remove tmp files
    os.remove( real_biom )
    os.remove( checked_biom )
