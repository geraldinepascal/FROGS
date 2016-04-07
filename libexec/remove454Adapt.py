#!/usr/bin/env python2.7
#
# Copyright (C) 2015 INRA
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
__version__ = '0.5.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import time
import argparse
import subprocess
from subprocess import Popen, PIPE

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
class TmpFiles:
    """
    @summary: Manager for temporary files.
    @note:
        tmpFiles = TmpFiles(out_dir)
        try:
            ...
            tmp_seq = tmpFiles.add( "toto.fasta" )
            ...
            tmp_log = tmpFiles.add( "log.txt" )
            ...
        finaly:
            tmpFiles.deleteAll()
    """
    def __init__(self, tmp_dir, prefix=None):
        """
        @param tmp_dir: [str] The temporary directory path.
        @param prefix: [str] The prefix added to each temporary file [default: <TIMESTAMP>_<PID>].
        """
        if prefix is None:
            prefix = str(time.time()) + "_" + str(os.getpid())
        self.files = list()
        self.tmp_dir = tmp_dir
        self.prefix = prefix

    def add(self, filename, prefix=None, dir=None):
        """
        @summary: Add a temporary file.
        @param filename: The filename without prefix.
        @param prefix: The prefix added [default: TmpFiles.prefix].
        @param dir: The directory path [default: TmpFiles.tmp_dir].
        @return: [str] The filepath.
        """
        # Default
        if prefix is None:
            prefix = self.prefix
        if dir is None:
            dir = self.tmp_dir
        # Process
        filepath = os.path.join(dir, prefix + "_" + filename)
        self.files.append(filepath)
        return filepath

    def delete(self, filepath):
        """
        @summary: Deletes the specified temporary file.
        @param filepath: [str] The file path to delete.
        """
        self.files.remove(filepath)
        if os.path.exists(filepath): os.remove(filepath)

    def deleteAll(self):
        """
        @summary: Deletes all temporary files.
        """
        all_tmp_files = [tmp_file for tmp_file in self.files]
        for tmp_file in all_tmp_files:
            self.delete(tmp_file)

def get_cutadapt_version():
    """
    @summary: Return the cutadapt version.
    @return: [str] The cutadapt version.
    """
    version = None
    try:
        cmd = 'cutadapt --version'
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        version = stdout.strip()
    except:
        version = "unknown"
    return version

def rvc( sequence ):
    """
    @summary: Return the reverse complement of the sequence.
    @param sequence: [str] The sequence to process.
    @return: [str] The reverse complement.
    """
    complement_rules = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N',
                        'a':'t','t':'a','g':'c','c':'g','u':'a','n':'n'}

    return "".join([complement_rules[base] for base in sequence[::-1]])


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Removes PCR primers on 454 amplicons.'
    )
    parser.add_argument( '-f', '--five-prim-primer', type=str, required=True, help="The 5' primer sequence (wildcards are accepted)." )
    parser.add_argument( '-t', '--three-prim-primer', type=str, required=True, help="The 3' primer sequence (wildcards are accepted)." )
    parser.add_argument( '-m', '--min-length', type=int, required=True, help="The amplicon minimal expected length (with primers)." )
    parser.add_argument( '-e', '--error-rate', type=float, default=0.1, help="Maximum allowed error rate (no. of errors divided by the length of the matching region)." )
    parser.add_argument( '-n', '--non-overlap', type=int, default=1, help="Maximum allowed error rate (no. of errors divided by the length of the matching region)." )
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ + " [cutadapt " + get_cutadapt_version() + "]" )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-i', '--input', required=True, help='The amplicons sequences (format: FASTQ).' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-o', '--output', required=True, help='The trimmed sequences.')
    args = parser.parse_args()
    if args.min_length < (len(args.five_prim_primer) + len(args.five_prim_primer)):
        raise argparse.ArgumentTypeError( "The minimum length of the amplicon (--min-length) must be superior to the size of the two primers." )

    # Process
    tmpFiles = TmpFiles( os.path.split(args.output)[0] )

    try:
        # Temporary files for user oriented reads
        tmp_5prim_trimmed = tmpFiles.add( "cutadapt_5prim_trimmed.fastq.gz" )
        tmp_5prim_untrimmed = tmpFiles.add( "cutadapt_5prim_not_trimmed.fastq" )
        tmp_5prim_too_short = tmpFiles.add( "cutadapt_5prim_too_short.fastq.gz" )
        tmp_5prim_log = tmpFiles.add( "cutadapt_5prim_log.txt" )
        tmp_3prim_trimmed = tmpFiles.add( "cutadapt_3prim_trimmed.fastq.gz" )
        tmp_3prim_log = tmpFiles.add( "cutadapt_3prim_log.txt" )
        # Temporary files for RVC oriented reads
        tmp_rvc_fastq = tmpFiles.add( "rvc.fastq.gz" )
        tmp_rvc_5prim_trimmed = tmpFiles.add( "cutadapt_rvc_5prim_trimmed.fastq.gz" )
        tmp_rvc_5prim_log = tmpFiles.add( "cutadapt_rvc_5prim_log.txt" )
        tmp_rvc_3prim_trimmed = tmpFiles.add( "cutadapt_rvc_3prim_trimmed.fastq.gz" )
        tmp_rvc_3prim_log = tmpFiles.add( "cutadapt_rvc_3prim_log.txt" )

        # Remove adapters
        subprocess.check_call( "cutadapt -g " + args.five_prim_primer +
                               " --error-rate " + str(args.error_rate) +
                               " --minimum-length " + str(args.min_length - len(args.five_prim_primer)) +
                               " --match-read-wildcards" +
                               " --overlap " + str(len(args.five_prim_primer) - args.non_overlap) + 
                               " --untrimmed-output " + tmp_5prim_untrimmed + " --too-short-output " + tmp_5prim_too_short + " -o " + tmp_5prim_trimmed +
                               " " + args.input +
                               " > " + tmp_5prim_log, shell=True) # Try to remove 5prim primer with sequences supposed to be forward
        subprocess.check_call( "cutadapt -a " + args.three_prim_primer +
                               " --error-rate " + str(args.error_rate) +
                               " --overlap " + str(len(args.three_prim_primer) - args.non_overlap) +
                               " --discard-untrimmed" +
                               " --match-read-wildcards" +
                               " -o " + tmp_3prim_trimmed +
                               " " + tmp_5prim_trimmed +
                               " > " + tmp_3prim_log, shell=True ) # For forward sequences (correctly trimmed by first cutadapt) remove 3prim primer

        # Collect supposed RVC reads (untrimmed and uncorrect trim)
        to_rvc = dict()
        for excluded in [tmp_5prim_untrimmed, tmp_5prim_too_short]:
            fh_in_excluded = FastqIO( excluded )
            for record in fh_in_excluded:
                to_rvc[record.id] = 1
            fh_in_excluded.close()

        fh_rvc_out = FastqIO( tmp_rvc_fastq, "w" )
        fh_initial_in = FastqIO( args.input )
        for record in fh_initial_in:
            if to_rvc.has_key( record.id ):
                record.string = rvc(record.string)
                fh_rvc_out.write( record )
        fh_initial_in.close()
        fh_rvc_out.close()

        # Remove adapters on RVC
        subprocess.check_call( "cutadapt -g " + args.five_prim_primer +
                               " --error-rate " + str(args.error_rate) +
                               " --discard-untrimmed" +
                               " --match-read-wildcards" +
                               " --overlap " + str(len(args.five_prim_primer) - args.non_overlap) +
                               "  -o " + tmp_rvc_5prim_trimmed +
                               " " + tmp_rvc_fastq +
                               " > " + tmp_rvc_5prim_log, shell=True )
        subprocess.check_call( "cutadapt -a " + args.three_prim_primer +
                               " --error-rate " + str(args.error_rate) +
                               " --discard-untrimmed" +
                               " --match-read-wildcards" +
                               " --overlap " + str(len(args.three_prim_primer) - args.non_overlap) +
                               " -o " + tmp_rvc_3prim_trimmed +
                               " " + tmp_rvc_5prim_trimmed +
                               " > " + tmp_rvc_3prim_log, shell=True )

        # Merge
        fh_final_out = FastqIO( args.output, "w" )
        for trimmed in [tmp_3prim_trimmed, tmp_rvc_3prim_trimmed]:
            fh_in_trimmed = FastqIO( trimmed )
            for record in fh_in_trimmed:
                fh_final_out.write( record )
            fh_in_trimmed.close()
        fh_rvc_out.close()
    finally:
        if not args.debug:
            tmpFiles.deleteAll()
