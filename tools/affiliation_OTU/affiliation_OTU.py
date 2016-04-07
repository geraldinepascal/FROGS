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

__author__ = 'Maria Bernard INRA - SIGENAE AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.10.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import json
import gzip
import argparse
import threading
import multiprocessing
import subprocess
from subprocess import Popen, PIPE

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
from frogsBiom import BiomIO
from frogsSequenceIO import *


###################################################################################################################
###                                 OTU AFFILIATION CLASSES                                                     ###
###################################################################################################################
class Blast(Cmd):
    def __init__(self, ref_fasta, query_fasta, output_blast, nb_cpus):
        """
        @param ref_fasta: [str] Path to the reference fasta file (blast indexed).
        @param query_fasta: [str] Path to the query fasta file to submit to blast
        @param output_blast: [str] Path to blast results.
        @param nb_cpus: [int] Number of usable CPUs.
        """
        Cmd.__init__( self,
                      "blastn",
                      "blast taxonomic affiliation",
                      "-num_threads " + str(nb_cpus) + " -task megablast -word_size 38 -max_target_seqs 500 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -query "+ query_fasta +" -out "+ output_blast +" -db " + ref_fasta,
                      "-version")

        self.output = output_blast

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').split()[1].strip()


class RDPAffiliation(Cmd):
    def __init__(self, ref, query_fasta, output, memory):
        """
        @param ref: [str] Path to the reference fasta file (rdp formated). This will automatically search the [ref].properties file
        @param query_fasta: [str] Path to the query fasta file to submit to RDP
        @param output: [str] Path to rdp results.
        @param memory: [int] Memory used by JVM (in Giga Bytes).
        """
        Cmd.__init__( self,
                      which("classifier.jar"),
                      "rdp taxonomic affiliation" ,
                      "taskset -c " + self._get_cpu_id() + " java -Xmx" + str(memory) + "g -jar ##PROGRAM## classify -c 0.0 -t " + ref + ".properties -o " + output + " " + query_fasta,
                      None)
        self.output = output

    def _get_cpu_id(self):
        python_pid = str(os.getpid())
        python_cpuid = None
        ps_stdout = subprocess.check_output('ps -o pid,cpuid', shell=True)
        for line in ps_stdout.strip().split("\n"):
            pid, cpuid = line.split()
            if pid == python_pid:
                python_cpuid = cpuid
        if python_cpuid is None:
            raise Exception( "CPUID cannot be retrieved" )
        return str(python_cpuid)


class AddAffiliation2Biom(Cmd):
    """
    @summary: Add Blast and/or RDP affiliation to biom
    @note: blast must be launched with -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'
    """
    def __init__(self, ref, blast_list, rdp_list, in_biom, out_biom):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_tsv: [str] Path to output TSV file.
        """
        argument="-f " + ref + " -i " + in_biom + " -o " + out_biom
        if blast_list is not None:
            argument += " -b " + " ".join(blast_list)
        if rdp_list is not None:
            argument += " -r " + " ".join(rdp_list)
        Cmd.__init__( self,
                      'addAffiliation2biom.py',
                      'Add Blast and/or RDP affiliation to biom',
                      argument,
                      '--version' )


###################################################################################################################
###                                 OTU AFFILIATION FUNCTIONS                                                   ###
###################################################################################################################
def get_fasta_nb_seq( fasta_file ):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(fasta_file):
        FH_input = open( fasta_file )
    else:
        FH_input = gzip.open( fasta_file )
    nb_seq = 0
    for line in FH_input:
        if line.startswith(">"):
            nb_seq += 1
    FH_input.close()
    return nb_seq

def split_fasta(fasta_file, tmp_files_manager, nb_file, out_list, log_file):
    """
    @summary: split fasta in nb_file and returne outfile list
    @param fasta_file: [str] Path to the fasta file to process.
    @param tmp_files_manager: [TmpFiles] The temporary file manager.
    @param nb_file : [int] The number of file to generate.
    @param out_list : [list] List of output fasta file
    @param log_file : [srt] path to logfile
    """
    out_files = list()
    record_iter = FastaIO(fasta_file)
    for idx, record in enumerate(record_iter):
        out_file_idx = idx % nb_file
        if len(out_files) == 0 or not out_file_idx < len(out_files):
            new_out_file = tmp_files_manager.add( os.path.basename(fasta_file) + "_" + str(out_file_idx) )
            out_files.append({ 'file_path': new_out_file,
                               'file_handle': FastaIO(new_out_file, "w"),
                               'nb_seq': 0
            })
            out_list.append( new_out_file )
        out_files[out_file_idx]['nb_seq'] += 1
        out_files[out_file_idx]['file_handle'].write( record )
    for out_file in out_files:
        out_file['file_handle'].close()

    # Log
    FH_log = Logger(log_file)
    FH_log.write("# split " + fasta_file + " in " + str(nb_file) + " fasta files\n")
    FH_log.write("Results\n")
    for out_file in out_files:
        FH_log.write( "\tWrote " + str(out_file['nb_seq']) + " records to " + out_file['file_path'] + "\n" )
    FH_log.close()

def process_rdp(input, output, log_file, reference, memory):
    """
    @summary: Launches RDP.
    @param input: [str] Path to the query fasta file to submit to RDP.
    @param output: [str] Path to rdp results.
    @param log_file: [str] Path to rdp log.
    @param reference: [str] Path to the reference fasta file (rdp formated).
    @param memory: [int] Memory used by RDP.
    """
    rdp_cmd = RDPAffiliation(reference, input, output, memory)
    rdp_cmd.submit(log_file)

def summarise_results( summary_file, biom_file, taxonomy_ranks ):
    """
    @summary: Writes one summary of results from several logs.
    @param summary_file: [str] The path to output file.
    @param biom_file: [str] The path to the BIOM file.
    @param taxonomy_ranks: [list] The ordered ranks levels present in the reference databank.
    """
    # Get data
    global_results, samples_results = get_results( biom_file )

    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "affiliation_OTU_tpl.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###GLOBAL_DATA###" in line:
            line = line.replace( "###GLOBAL_DATA###", json.dumps(global_results) )
        elif "###SAMPLES_DATA###" in line:
            line = line.replace( "###SAMPLES_DATA###", json.dumps(samples_results) )
        elif "###TAXONOMY_RANKS###" in line:
            line = line.replace( "###TAXONOMY_RANKS###", json.dumps(taxonomy_ranks) )
        FH_summary_out.write( line )
    FH_summary_out.close()
    FH_summary_tpl.close()

def get_results( biom_file ):
    """
    @summary: Returns the results of the affiliation.
    @param biom_file: [str] Path to a BIOM file after affiliation.
    @return: [dict] The global results and the sample results.
    """
    global_results = {
        "nb_clstr": 0,
        "nb_seq": 0,
        "nb_clstr_with_affi": 0,
        "nb_seq_with_affi": 0,
        "nb_clstr_ambiguous": list(),
        "nb_seq_ambiguous": list(),
    }
    samples_results = dict()

    biom = BiomIO.from_json( biom_file )
    for cluster in biom.get_observations():
        nb_seq = biom.get_observation_count( cluster["id"] )
        global_results["nb_clstr"] += 1
        global_results["nb_seq"] += nb_seq
        if cluster["metadata"]["blast_taxonomy"] is not None:
            global_results["nb_clstr_with_affi"] += 1
            global_results["nb_seq_with_affi"] += nb_seq
            for depth, taxon in enumerate(cluster["metadata"]["blast_taxonomy"]):
                if len(global_results["nb_clstr_ambiguous"]) < (depth + 1):
                    global_results["nb_clstr_ambiguous"].append( 0 )
                    global_results["nb_seq_ambiguous"].append( 0 )
                if taxon == "Multi-affiliation":
                    global_results["nb_clstr_ambiguous"][depth] += 1
                    global_results["nb_seq_ambiguous"][depth] += nb_seq
        # Samples results
        for sample in biom.get_samples_by_observation( cluster["id"] ):
            sample_name = sample["id"]
            if not samples_results.has_key( sample_name ):
                samples_results[sample_name] = {
                    "nb_clstr": 0,
                    "nb_seq": 0,
                    "nb_clstr_with_affi": 0,
                    "nb_seq_with_affi": 0
                }
            count = biom.get_count(cluster["id"], sample_name)
            if count > 0:
                samples_results[sample_name]["nb_clstr"] += 1
                samples_results[sample_name]["nb_seq"] += count
                if cluster["metadata"]["blast_taxonomy"] is not None:
                    samples_results[sample_name]["nb_clstr_with_affi"] += 1
                    samples_results[sample_name]["nb_seq_with_affi"] += count
    return global_results, samples_results


def log_append_files(log_file, appended_files):
    """
    @summary: Append content of several log files in one log file.
    @param log_file: [str] The log file where contents of others are appended.
    @param appended_files: [list] List of log files to append.
    """
    FH_log = Logger(log_file)
    FH_log.write("\n")
    for current_file in appended_files:
        FH_input = open(current_file)
        for line in FH_input:
            FH_log.write(line)
        FH_input.close()
        FH_log.write("\n")
    FH_log.write("\n")
    FH_log.close()

def parallel_submission( function, inputs, outputs, logs, cpu_used, reference, memory):
    processes = [{'process':None, 'inputs':None, 'outputs':None, 'log_files':None} for idx in range(cpu_used)]
    # Launch processes
    for idx in range(len(inputs)):
        process_idx = idx % cpu_used
        processes[process_idx]['inputs'] = inputs[idx]
        processes[process_idx]['outputs'] = outputs[idx]
        processes[process_idx]['log_files'] = logs[idx]

    for current_process in processes:
        if idx == 0:  # First process is threaded with parent job
            current_process['process'] = threading.Thread(target=function,
                                                          args=(current_process['inputs'], current_process['outputs'], current_process['log_files'], reference, memory))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                                 args=(current_process['inputs'], current_process['outputs'], current_process['log_files'], reference, memory))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise Exception("Error in sub-process execution.")


###################################################################################################################
###                                              MAIN                                                           ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Taxonomic affiliation of each OTU's seed by RDPtools and BLAST.")
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used.(default 1)")
    parser.add_argument( '-m', '--java-mem', type=int, default=2, help="Java memory allocation in Go.(default 2)")
    parser.add_argument( '-t', '--taxonomy-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels present in the reference databank.' )
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-r', '--reference', required=True, help='Preformated reference file.')
    group_input.add_argument('-b', '--input-biom', required=True, help='Abundance table from the clusterisation program (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help="Fasta file of OTU's seed (format: fasta).")
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-biom', default='affiliation.biom', help='File which add affiliation annotations from blast and RDPtools to the abundance table.')
    group_output.add_argument('-s', '--summary', default='summary.html', help='Report of the results (format: HTML).')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )
    fasta_rdp_list = []
    rdp_out_list = []
    log_rdp_list = []
    fasta_blast_list = []
    blast_out_list = []
    log_blast_list = []

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
        nb_seq = get_fasta_nb_seq(args.input_fasta)
        Logger.static_write(args.log_file, "Nb seq : " + str(nb_seq) + "\n\n")

        if args.nb_cpus == 1 or nb_seq < 100:
            # RDP
            rdp_out_list.append( tmpFiles.add(os.path.basename(args.input_fasta) + ".rdp") )
            process_rdp( args.input_fasta, rdp_out_list[0], args.log_file, args.reference, args.java_mem)
            # Blast
            blast_out_list.append( tmpFiles.add(os.path.basename(args.input_fasta) + ".blast") )
            Blast(args.reference, args.input_fasta, blast_out_list[0], 1).submit(args.log_file)
        else:
            # RDP
            split_fasta(args.input_fasta, tmpFiles, max(1, int(args.nb_cpus/3)), fasta_rdp_list, args.log_file)
            rdp_out_list = [tmpFiles.add(os.path.basename(current_fasta) + ".rdp") for current_fasta in fasta_rdp_list]
            log_rdp_list = [tmpFiles.add(os.path.basename(current_fasta) + "_rdp.log") for current_fasta in fasta_rdp_list]
            parallel_submission( process_rdp, fasta_rdp_list, rdp_out_list, log_rdp_list, len(fasta_rdp_list), args.reference, args.java_mem )
            # Blast
            blast_out_list.append( tmpFiles.add(os.path.basename(args.input_fasta) + ".blast") )
            log_blast_list.append( tmpFiles.add(os.path.basename(args.input_fasta) + "_blast.log") )
            Blast(args.reference, args.input_fasta, blast_out_list[0], args.nb_cpus).submit(log_blast_list[0])
            # Logs
            log_append_files(args.log_file, log_rdp_list + log_blast_list)

        # Convert to output file
        AddAffiliation2Biom( args.reference, blast_out_list, rdp_out_list, args.input_biom, args.output_biom ).submit( args.log_file )
        summarise_results( args.summary, args.output_biom, args.taxonomy_ranks )

    finally:
        if not args.debug:
            tmpFiles.deleteAll()