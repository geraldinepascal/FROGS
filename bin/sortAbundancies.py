#!/usr/bin/env python2.7
#
# Copyright (C) 2014 INRA
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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.4.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import time
import gzip
import argparse
from sequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def sort_by_abundance( input_fasta, output_fasta, tmp_folder, ranges, old_size_separator, new_size_separator="_", sort_type="desc" ):
    prefix = str(time.time()) + "_" + str(os.getpid())
    abundances_files = dict()

    try:
        # Split by abundance
        FH_input = FastaIO( input_fasta )
        for record in FH_input:
            abundance = 1
            abundance_max = 1
            id_without_abundance = record.id
            if old_size_separator != None and old_size_separator in record.id:
                abundance = int(record.id.rsplit(old_size_separator)[-1])
                idx = 0
                while abundance > ranges[idx]:
                    idx += 1
                abundance_max = ranges[idx]
                id_without_abundance = record.id.rsplit(old_size_separator, 1)[0]
            if not abundances_files.has_key( abundance_max ):
                abundances_files[abundance_max] = dict()
                abundances_files[abundance_max]["path"] = os.path.join(tmp_folder, prefix + '_tmp_abundance_' + str(abundance_max) +  '.fasta' )
                abundances_files[abundance_max]["filehandle"] = FastaIO( abundances_files[abundance_max]["path"], "w" )
            record.id = id_without_abundance + new_size_separator + str(abundance)
            abundances_files[abundance_max]["filehandle"].write( record )
        for abundance in abundances_files.keys():
            abundances_files[abundance]["filehandle"].close()
        # Merge
        output_fh = FastaIO( output_fasta, "w" )
        reverse_order = True if sort_type == "desc" else False
        abundances = sorted( abundances_files.keys(), key=int, reverse=reverse_order )
        for idx, current_abundance in enumerate(abundances):
            if (idx != 0 and abs(current_abundance - abundances[idx-1]) == 1) or (idx == 0 and current_abundance == 1 and not reverse_order): # File with only one abundance
                fh = FastaIO( abundances_files[current_abundance]["path"] )
                for current_record in fh:
                    output_fh.write( current_record )
                fh.close()
            else: # File with different abundances
                records_by_abundance = dict()
                fh = FastaIO( abundances_files[current_abundance]["path"] )
                for record in fh:
                    abundance = int( record.id.rsplit(new_size_separator, 1)[-1] )
                    if not records_by_abundance.has_key(abundance):
                        records_by_abundance[abundance] = list()
                    records_by_abundance[abundance].append(record)
                fh.close()
                for current_record_abundance in sorted( records_by_abundance.keys(), key=int, reverse=reverse_order ):
                    for current_record in records_by_abundance[current_record_abundance]:
                        output_fh.write( current_record )
                del records_by_abundance
        output_fh.close()
    finally:
        #Remove tmp
        for size in abundances_files.keys():
            current_file = abundances_files[size]["path"]
            if os.path.exists(current_file) : os.remove( current_file )

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Sort fasta y abundancies.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-s', '--size-separator', default=None, help='Each sequence in in_fasta is see as a pre-cluster. The number of sequences represented by the pre-cluster is stored in sequence ID. Sequence ID format : "<REAL_ID><size_separator><NB_SEQ>". If this size separator is missing in ID, the number of sequences represented is 1.' )
    parser.add_argument( '-i', '--input-file', required=True, help='Fasta file to sort.' )
    parser.add_argument( '-o', '--output-file', required=True, help='Sorted fasta file.' )
    parser.add_argument( '-r', '--ranges', type=int, nargs='*', default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 100, 200, 300, 400, 500, 600, 700, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 100000, 1000000, 9999999999999], help='Ranges used to decrease the memory comsumption. Each values represents one file which contains the sequences with an abundance between previous value and next values. Each file is processed separately.' )
    args = parser.parse_args()

    # Process
    tmp_folder = os.path.split(args.output_file)[0]
    sort_by_abundance( args.input_file, args.output_file, tmp_folder, sorted(args.ranges), args.size_separator )
