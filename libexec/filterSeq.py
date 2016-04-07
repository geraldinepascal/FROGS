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
import gzip
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from frogsSequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
class QualWindowParameter(argparse.Action):
    """
    @summary : Argparse parameter for qual_window.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        valid_subparameters = ['win_size', 'offset', 'threshold']
        # Set parser
        param = getattr(namespace, self.dest)
        if param is None:
            param = {
                'win_size': 20,
                'offset': 33
            }
        # Retrieve params
        for option in values:
            key, val = option.split(":")
            if key not in valid_subparameters:
                raise argparse.ArgumentTypeError( "For --qual-window only the following sub-parameters are authorized: " + ", ".join(valid_subparameters) )
            param[key] = int(val)
        setattr(namespace, self.dest, param)

def filter_seq( input_file, output_file, log_file, min_length=None, max_length=None, max_N=None, max_homopolymer=None, qual_window=None, force_fasta=False ):
    """
    @summary: Filters sequences on length, number of N and number of homopolymer.
    @param input_file: [str] Path to the processed sequence file.
    @param output_file: [str] Path to the output file.
    @param log_file: [str] Path to the log file. It contains the numeric results of each filters.
    @param min_length: [str] The minimum length to keep a sequence. The value None cancels the filter on this criterion.
    @param max_length: [str] The maximum length to keep a sequence. The value None cancels the filter on this criterion.
    @param max_N: [str] The maximum number of N to keep a sequence. The value None cancels the filter on this criterion.
    @param max_homopoly: [str] The maximum number of homopolymer to keep a sequence. The value None cancels the filter on this criterion.
    @param qual_window: [dict] The minimal distance between two poor quality.
    @param force_fasta: [bool] With True the output will be a fasta. Otherwise the output format is the same of the input.
    """
    def is_true( *val ):
        return True

    def check_min_length( length ):
        return length >= int(min_length)

    def check_max_length( length ):
        return length <= int(max_length)

    def check_length( length ):
        return( length >= int(min_length) and length <= int(max_length) )

    def check_N_number( sequence ):
        return( (sequence.count('N') + sequence.count('n')) <= int(max_N) )

    def check_homopolymer( sequence ):
        previous = ''
        homopolymer_length = 0
        sequence_length = len( sequence )
        idx = 0
        while homopolymer_length < int(max_homopolymer) + 1 and idx < sequence_length:
            if sequence[idx] == previous:
                homopolymer_length += 1
            else:
                homopolymer_length = 0
            previous = sequence[idx]
            idx += 1
        return homopolymer_length < int(max_homopolymer)

    def check_quality( quality ):
        is_keep = True
        previous_idx = -qual_window['win_size'] - 1
        for idx in range(len(quality)):
            if (ord(quality[idx]) - qual_window['offset']) < qual_window['threshold']:
                if (idx - previous_idx - 1) <= qual_window['win_size']:
                    is_keep = False
                previous_idx = idx
        return is_keep

    # Load suitable function to check length
    length_is_ok = is_true # no check
    if min_length is not None and max_length is not None: # check min and max length
        length_is_ok = check_length
    elif min_length is not None and max_length is None: # check only min length
        length_is_ok = check_min_length
    elif min_length is None and max_length is not None: # check only max length
        length_is_ok = check_max_length

    # Load suitable function to check the number of N
    N_number_is_ok = is_true # no check
    if max_N is not None:
        N_number_is_ok = check_N_number

    # Load suitable function to check the number homopolymer
    homopolymer_is_ok = is_true # no check
    if max_homopolymer is not None:
        homopolymer_is_ok = check_homopolymer

    # Load suitable function to check the quality
    quality_is_ok = is_true # no check
    if qual_window is not None:
        quality_is_ok = check_quality

    # Process
    #~ fh_in = SequenceFileReader.factory(input_file)
    fh_in = FastqIO(input_file)
    if force_fasta:
        fh_out = FastaIO(output_file, "w")
    elif issubclass(fh_in.__class__, FastqIO):
        fh_out = FastqIO(output_file, "w")
    else:
        fh_out = FastaIO(output_file, "w")
    nb_seq = 0
    filter_on_length = 0
    filter_on_N = 0
    filter_on_homopoly = 0
    filter_on_quality = 0
    for seq_record in fh_in:
        nb_seq += 1
        if not length_is_ok( len(seq_record.string) ):
            filter_on_length += 1
        elif not N_number_is_ok( seq_record.string ):
            filter_on_N += 1
        elif not homopolymer_is_ok( seq_record.string ):
            filter_on_homopoly += 1
        elif not quality_is_ok( seq_record.quality ):
            filter_on_quality += 1
        else:
            fh_out.write( seq_record )

    # Write log
    log_fh = open( log_file, "w" )
    log_fh.write( "Nb seq processed : " + str(nb_seq) + "\n" )
    if not(min_length is None and max_length is None):
        log_fh.write( "Nb seq filtered on length : " + str(filter_on_length) + "\n" )
    if not max_N is None:
        log_fh.write( "Nb seq filtered on N : " + str(filter_on_N) + "\n" )
    if not max_homopolymer is None:
        log_fh.write( "Nb seq filtered on homopolymer : " + str(filter_on_homopoly) + "\n" )
    if not qual_window is None:
        log_fh.write( "Nb seq filtered on quality : " + str(filter_on_quality) + "\n" )
    log_fh.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Filter sequences on length, number of N, number of homopolymer and distance between positions with a poor quality.' )
    parser.add_argument( '--min-length', default=None, type=int, help='The minimum length to keep a sequence.' )
    parser.add_argument( '--max-length', default=None, type=int, help='The maximum length to keep a sequence.' )
    parser.add_argument( '--max-N', default=None, type=int, help='The maximum number of N to keep a sequence.' )
    parser.add_argument( '--max-homopolymer', default=None, type=int, help='The maximum number of homopolymer to keep a sequence.' )
    parser.add_argument( '--qual-window', default=None, nargs='+', action=QualWindowParameter, metavar=("threshold:MIN_QUAL [win_size:WINDOW_SIZE] [offset:QUAl_OFFSET]", ""), help="The minimal distance ('win_size') between two poor quality ('threshold')." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file', required=True, help='The sequence file to process (format : FASTQ or FASTA).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file', required=True, help='The filtered file.')
    group_output.add_argument( '-l', '--log-file', required=True, help='The log file.')
    group_output.add_argument( '-c', '--compress', default=False, action='store_true', help='Compress the output file (algorithm : gzip). If necessary the script add ".gz" suffix to output name.' )
    group_output.add_argument( '-f', '--force-fasta', default=False, action='store_true', help='The output will be a fasta. Otherwise the output format is the same of the input.' )
    args = parser.parse_args()

    output_file = args.output_file
    if args.compress:
        if not args.output_file.endswith('gz'):
            output_file = args.output_file + '.gz'
    filter_seq( args.input_file, output_file, args.log_file, args.min_length, args.max_length, args.max_N, args.max_homopolymer, args.qual_window, args.force_fasta )