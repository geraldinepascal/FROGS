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
__version__ = '0.2.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'beta'

import os
import sys
import time
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR)
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR

from subprocess import Popen, PIPE
from frogsBiom import BiomIO
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

def submit_cmd( cmd, stdout_path, stderr_path ):
    """
    @summary: Submits the command and checks its exit status.
    @param cmd: [list] The command.
    @param stdout_path: [str] Filepath to the stdout.
    @param stderr_path: [str] Filepath to the stderr.
    """
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # write down the stdout
    stdoh = open(stdout_path, "w")
    stdoh.write(stdout)
    stdoh.close()

    # write down the stderr
    stdeh = open(stderr_path, "w")
    stdeh.write(stderr)
    stdeh.close()

    # check error status
    if p.returncode != 0:
        stdeh = open(stderr_path)
        error_msg = "".join( map(str, stdeh.readlines()) )
        stdeh.close()
        raise StandardError( error_msg )

def get_conta( min_identity, min_coverage, blast_file ):
    """
    @summary: Returns the ID of contaminated sequences.
    @param min_identity: [float] The minimum identity between sequence and contaminant for a contaminated sequence.
    @param min_coverage: [float] The minimum query coverage between sequence and contaminant for a contaminated sequence.
    @param blast_file: [str] Path to the alignment between contaminants and sequences.
    @return: [list] the contaminated_ids (dict) and the alignment_heatmap (list of coverage by identity).
    """
    contaminated_ids = dict()
    alignment_heatmap = [[0 for coverage in range(101)] for identity in range(101)]

    # Retrieve contaminated sequences
    current_hit = None
    FH_blast = open(blast_file)
    for line in open(blast_file):
        parts = line.strip().split()
        query_id = parts[0].split(";")[0]
        hit_id = parts[1]
        score = parts[11]
        perc_identity = float(parts[2])
        perc_coverage = ((int(parts[7]) - int(parts[6]) + 1)/float(parts[12]))*100
        if current_hit is not None and query_id + hit_id != current_hit['query_id'] + current_hit['hit_id']:
            # Store alignment info
            alignment_heatmap[int(current_hit['perc_identity']+0.5)][int(current_hit['perc_coverage']+0.5)] += 1
            # Filter alignment
            if current_hit['perc_identity'] >= float(min_identity)*100 and current_hit['perc_coverage'] >= float(min_coverage)*100:
                contaminated_ids[current_hit['query_id']] = True
            # New hit
            current_hit = None
        if current_hit is None or (score > current_hit['best_score']):
            current_hit = {
                'hit_id': hit_id,
                'query_id': query_id,
                'best_score': score,
                'perc_identity': perc_identity,
                'perc_coverage': perc_coverage
            }
    FH_blast.close()
    if current_hit is not None:
        # Store alignment info
        alignment_heatmap[int(current_hit['perc_identity']+0.5)][int(current_hit['perc_coverage']+0.5)] += 1
        # Filter alignment
        if current_hit['perc_identity'] >= float(min_identity)*100 and current_hit['perc_coverage'] >= float(min_coverage)*100:
            contaminated_ids[current_hit['query_id']] = True

    return contaminated_ids, alignment_heatmap

def write_log( log_file, processed_count, removed_count, alignment_heatmap ):
    """
    @summary: Write log.
    @param log_file: [str] Path to the log file.
    @param processed_count: [int] The number of processed sequences.
    @param removed_count: [int] The number of removed sequences.
    @param alignment_heatmap: [list] The alignment heatmap (coverage by identity).
    """
    # Write logs
    FH_log = open(log_file, "w")
    FH_log.write("##Global info\n")
    FH_log.write("#Processed :\t" + str(processed_count) + "\n")
    FH_log.write("#Contaminated :\t" + str(removed_count) + "\n")
    FH_log.write("#Clean :\t" + str(processed_count - removed_count) + "\n")
    FH_log.write("\n")
    FH_log.write("##Alignment info\n")
    header = "#Identity|Coverage\t" + "\t".join([str(coverage) for coverage in range(101)])
    FH_log.write(header + "\n")
    for identity, coverages in enumerate(alignment_heatmap):
        line = str(identity) + "\t" + "\t".join(map(str, coverages))
        FH_log.write(line + "\n")
    FH_log.close()

def filter_biom( removed_observations, in_biom, out_biom ):
    """
    @summary: Removed the specified observations from BIOM.
    @param removed_observations: [dict] Each key is an observation name.
    @param in_biom: [str]: Path to the processed BIOM file.
    @param out_biom: [str]: Path to the cleaned BIOM file.
    """
    biom = BiomIO.from_json(in_biom)
    biom.remove_observations(removed_observations)
    BiomIO.write(out_biom, biom)

def split_fasta( contaminated, in_fasta, cleaned_fasta, removed_fasta ):
    """
    @summary: Splits contaminated sequences and other sequences in two separated files.
    @param input_fasta: [str] Path to the processed fasta file.
    @param cleaned_fasta: [str] Path to the output file without contaminated sequences.
    @param removed_fasta: [str] Path to the output file with only contaminated sequences.
    """
    nb_seq = 0

    FH_total = FastaIO(in_fasta)
    FH_cleaned = FastaIO(cleaned_fasta, "w")
    FH_removed = FastaIO(removed_fasta, "w")
    for record in FH_total:
        nb_seq += 1
        if contaminated.has_key(record.id):
            FH_removed.write( record )
        else:
            FH_cleaned.write( record )
    FH_total.close()
    FH_cleaned.close()
    FH_removed.close()

    return nb_seq


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Uses similarity with a contaminant databank to split contaminated sequences in first file and others sequences in an other file.")
    parser.add_argument('-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used.(default 1)")
    parser.add_argument('-w', '--word-size', type=int, default=40, help="Word size for blast wordfinder algorithm (length of best perfect match).")
    parser.add_argument('--min-identity', type=float, default=0.8, help="Minimum identity between query and databank sequence to tag the query as contaminant.")
    parser.add_argument('--min-coverage', type=float, default=0.8, help="Minimum coverage between query and databank sequence to tag the query as contaminant.")
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument( '-f', '--input-fasta', required=True, help='The sequences to check (format: FASTA).' )
    group_input.add_argument( '-c', '--contaminant-db', required=True, help='The sequences of contaminants (format: FASTA with blast index).')
    group_input.add_argument( '-b', '--input-biom', default=None, help='The abundance file (format: BIOM).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument( '--clean-fasta', default='clean.fasta', help='The sequences without contaminants (format: FASTA).' )
    group_output.add_argument( '--clean-biom', default=None, help='The abundance without contaminants (format: BIOM).' )
    group_output.add_argument( '--conta-fasta', default=None, help='The checked sequences tagged as contaminants (format: FASTA).' )
    group_output.add_argument( '-l', '--log-file', default='log.txt', help='The log file.' )
    args = parser.parse_args()

    tmp_files = TmpFiles(os.path.split(args.clean_fasta)[0])
    try:
        blast_output = tmp_files.add("blast.tsv")
        blast_stdout = tmp_files.add("blast.stdout")
        blast_stderr = tmp_files.add("blast.stderr")
        contaminated_fasta = args.conta_fasta if args.conta_fasta is not None else tmp_files.add("contaminated.fasta")

        # Blast
        blast_cmd = ["blastn", "-num_threads", str(args.nb_cpus), "-word_size", str(args.word_size), "-max_target_seqs", "1",
                     "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
                     "-query", args.input_fasta, "-out", blast_output, "-db", args.contaminant_db]
        submit_cmd( blast_cmd, blast_stdout, blast_stderr )

        # Split contaminants
        contaminated_ids, alignment_heatmap = get_conta( args.min_identity, args.min_coverage, blast_output )
        nb_seq_ini = split_fasta( contaminated_ids, args.input_fasta, args.clean_fasta, contaminated_fasta )
        if args.input_biom is not None and args.clean_biom:
            filter_biom( contaminated_ids, args.input_biom, args.clean_biom )

        # Log
        write_log( args.log_file, nb_seq_ini, len(contaminated_ids), alignment_heatmap )
    finally:
        if not args.debug:
            tmp_files.deleteAll()