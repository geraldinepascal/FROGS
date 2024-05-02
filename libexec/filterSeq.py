#!/usr/bin/env python3
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
__version__ = '1.6.0'
__email__ = 'frogs-support@inrae.fr'
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
from frogsUtils import *


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
                raise_exception( argparse.ArgumentTypeError( "\n\n#ERROR : For --qual-window only the following sub-parameters are authorized: " + ", ".join(valid_subparameters) + '\n\n'))
            param[key] = int(val)
        setattr(namespace, self.dest, param)

def filter_seq( input_file1, output_file1, log_file, min_length=None, max_length=None, max_N=None, max_homopolymer=None, qual_window=None, tag=None, force_fasta=False, input_file2=None, output_file2=None):
    """
    @summary: Filters sequences on length, number of N and number of homopolymer.
    @param input_file1: [str] Path to the processed sequence file1.
    @param output_file1: [str] Path to the output file1.
    @param log_file: [str] Path to the log file. It contains the numeric results of each filters.
    @param min_length: [str] The minimum length to keep a sequence. The value None cancels the filter on this criterion.
    @param max_length: [str] The maximum length to keep a sequence. The value None cancels the filter on this criterion.
    @param max_N: [str] The maximum number of N to keep a sequence. The value None cancels the filter on this criterion.
    @param max_homopoly: [str] The maximum number of homopolymer to keep a sequence. The value None cancels the filter on this criterion.
    @param qual_window: [dict] The minimal distance between two poor quality.
    @param tag: [str] A particular tag included in the sequence to keep it. The value None cancels the filter on this criterion.
    @param force_fasta: [bool] With True the output will be a fasta. Otherwise the output format is the same of the input.
    @param input_file2: [str] Path to the processed sequence file2.
    @param output_file2: [str] Path to the output file2.
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

    def check_tag(sequence):
        return (tag in sequence)

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

    # Load suitable function to check the presence of tag
    tag_seq_is_ok = is_true # no check
    if tag is not None:
        tag_seq_is_ok = check_tag

    # Process
    fh_in = SequenceFileReader.factory(input_file1)
    if input_file2 is not None:
        fh_in2 = SequenceFileReader.factory(input_file2)

    # fh_in = FastqIO(input_file)
    if force_fasta:
        fh_out = FastaIO(output_file1, "wt")
        if input_file2 is not None:
            fh_out2 = FastaIO(output_file2, "wt")
    elif issubclass(fh_in.__class__, FastqIO):
        fh_out = FastqIO(output_file1, "wt")
        if input_file2 is not None:
            fh_out2 = FastqIO(output_file2, "wt")
    else:
        fh_out = FastaIO(output_file1, "wt")
        if input_file2 is not None:
            fh_out = FastaIO(output_file2, "wt")
    nb_seq = 0
    filter_on_length = 0
    filter_on_N = 0
    filter_on_homopoly = 0
    filter_on_quality = 0
    filter_on_tag = 0
    for seq_record in fh_in:
        if input_file2 is not None:
            seq_record2=fh_in2.next_seq()
        
        if ";size=" in seq_record.id :
            cur_nb_seq = int(seq_record.id.split(';size=')[1])
        else:
            cur_nb_seq = 1
        nb_seq += cur_nb_seq
        
        if input_file2 is None:
            if not tag_seq_is_ok(seq_record.string):
                filter_on_tag += cur_nb_seq
            elif not length_is_ok( len(seq_record.string) ):
                filter_on_length += cur_nb_seq
            elif not N_number_is_ok( seq_record.string ): 
                filter_on_N += cur_nb_seq
            elif not homopolymer_is_ok( seq_record.string ):
                filter_on_homopoly += cur_nb_seq
            elif not quality_is_ok( seq_record.quality ):
                filter_on_quality += cur_nb_seq
            else:
                fh_out.write( seq_record )
        else:
            print(len(seq_record.string))
            if not tag_seq_is_ok(seq_record.string):
                filter_on_tag += cur_nb_seq
            elif not N_number_is_ok( seq_record.string ) or not N_number_is_ok( seq_record2.string ):
                filter_on_N += cur_nb_seq
            elif not length_is_ok( len(seq_record.string) ) or not length_is_ok( len(seq_record2.string )):
                filter_on_length += cur_nb_seq
            else:
                fh_out.write( seq_record )
                fh_out2.write( seq_record2 )

    # Write log
    log_fh = open( log_file, "wt" )
    log_fh.write( "Nb seq processed : " + str(nb_seq) + "\n" )
    if not(min_length is None or min_length == 20 and max_length is None):
        log_fh.write( "Nb seq filtered on length : " + str(filter_on_length) + "\n" )
    if not max_N is None:
        log_fh.write( "Nb seq filtered on N : " + str(filter_on_N) + "\n" )
    if not max_homopolymer is None:
        log_fh.write( "Nb seq filtered on homopolymer : " + str(filter_on_homopoly) + "\n" )
    if not qual_window is None:
        log_fh.write( "Nb seq filtered on quality : " + str(filter_on_quality) + "\n" )
    if not tag is None:
        log_fh.write( "Nb seq filtered on absence of tag : " + str(filter_on_tag) + "\n" )
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
    parser.add_argument( '--tag', default=None, type=str, help='A particular tag included in the sequence to keep it' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '-i', '--input-file1', required=True, help='The (R1) sequence file to process (format: fastq or fasta).' )
    group_input.add_argument( '-j', '--input-file2', help='The (R2) sequence file to process (format: fastq or fasta).' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '-o', '--output-file1', required=True, help='The filtered sequence file1 (R1).')
    group_output.add_argument( '-p', '--output-file2', help='The filtered sequence file2 (R2).')
    group_output.add_argument( '-l', '--log-file', required=True, help='The log file.')
    group_output.add_argument( '-c', '--compress', default=False, action='store_true', help='Compress the output file (algorithm : gzip). If necessary the script add ".gz" suffix to output name.' )
    group_output.add_argument( '-f', '--force-fasta', default=False, action='store_true', help='The output will be a fasta. Otherwise the output format is the same of the input.' )
    args = parser.parse_args()

    if args.compress:
        if not args.output_file1.endswith('gz'):
            output_file1 = args.output_file1 + '.gz'
            output_file2 = args.output_file2 + '.gz'
    filter_seq( args.input_file1, args.output_file1, args.log_file, args.min_length, args.max_length, args.max_N, args.max_homopolymer, args.qual_window, args.tag, args.force_fasta, args.input_file2, args.output_file2)
