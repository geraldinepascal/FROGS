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

__author__ = 'Maria Bernard INRA - SIGENAE AND Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '3.2'
__email__ = 'frogs-support@inrae.fr'
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
else: os.environ['PYTHONPATH'] = LIB_DIR + os.pathsep + os.environ['PYTHONPATH']

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

        # self.output = output_blast

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout').split()[1].strip()

class SplitOnTag(Cmd):
    """
    @summary : split fasta sequence on tag
    """
    def __init__(self, combined_input , split_tag , out_split_1 , out_split_2 ):
        """
        @param combined_input  : [str] Path to combined sequence 
        @param split_tag     : [str] the sequence tag on which to split sequences
        @param out_split_1     : [str] Path to fasta split sequence output file 1
        @param out_split_2     : [str]  Path to fasta split sequence output file 2
        """
        Cmd.__init__( self,
                      'combine_and_split.py',
                      'Split on tag.',
                      ' --reads1 ' + combined_input + ' -s ' + split_tag + ' --split-output1 ' + out_split_1 + ' --split-output2 ' + out_split_2,
                      '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()

class Reduce_Ref(Cmd):
    def __init__(self,input_ref, query_blast_R1, query_blast_R2, output_fasta, log):
        
        Cmd.__init__( self,
                      "reduce_ref_for_needleall.py",
                      "Reduced reference fasta file based on R1 and R2 100 best score blast alignment",
                      " -r " + input_ref + " --query-blast-R1 " + query_blast_R1 + " --query-blast-R2 " + query_blast_R2 + " -f " + output_fasta + " -l " + log,
                      "-v")
        
        self.reduced_log=log

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').strip()
        
        
    def parser(self, log_file):

        FH_in = open(self.reduced_log)
        FH_log = Logger( log_file )
        for line in FH_in:
            FH_log.write( '\t'+line )
        FH_log.write('\n')
        FH_log.close()
        FH_in.close()


class Needleall(Cmd):
    def __init__(self,input_ref, input_query, output_sam, log):
    
        Cmd.__init__( self,
                  "needleall",
                  "Perform global alignment",
                  " -asequence " + input_ref + " -bsequence " + input_query + " -outfile " + output_sam + " -aformat3 sam -gapopen 10.0 -gapextend 0.5 -errfile " + log,
                  "--version")
                  
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').strip()

class NeedleallSam_to_tsv(Cmd):
    def __init__(self, input_sam, input_ref, output_tsv):
        """
        @param input_sam: [str] Path to NeedleAll Sam output file.
        @param input_ref: [str] Path to reference fasta file use for NeedleAll alignment.
        @param output_tsv: [int] Path to ouput Blast like tsv file.
        """
        Cmd.__init__( self,
                      "needleallSam_to_tsv.py",
                      "convert NeedleAll Sam output in blast like tsv output sorted by bitscore",
                      "-n " + input_sam + " -r " + input_ref + " -b " + output_tsv,
                      "-v")

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stderr').strip()

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
                      #~ "taskset -c " + self._get_cpu_id() + " java -Xmx" + str(memory) + "g -jar ##PROGRAM## classify -c 0.0 -t " + ref + ".properties -o " + output + " " + query_fasta,
                      "java -Xmx" + str(memory) + "g -jar ##PROGRAM## classify -c 0.0 -t " + ref + ".properties -o " + output + " " + query_fasta,
                      None)
        self.output = output

    # def _get_cpu_id(self):
    #     python_pid = str(os.getpid())
    #     python_cpuid = None
    #     ps_stdout = subprocess.check_output('ps -o pid,cpuid', shell=True)
    #     for line in ps_stdout.strip().split("\n"):
    #         pid, cpuid = line.split()
    #         if pid == python_pid:
    #             python_cpuid = cpuid
    #     if python_cpuid is None:
    #         raise Exception( "\nCPUID cannot be retrieved\n\n" )
    #     return str(python_cpuid)


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
        if len(blast_list) > 0 :
            argument += " -b " + " ".join(blast_list)
        if len(rdp_list) > 0 :
            argument += " -r " + " ".join(rdp_list)
        Cmd.__init__( self,
                      'addAffiliation2biom.py',
                      'Add Blast and/or RDP affiliation to biom',
                      argument,
                      '--version' )


###################################################################################################################
###                                 OTU AFFILIATION FUNCTIONS                                                   ###
###################################################################################################################
def extract_FROGS_combined(input_fasta, fasta_full_length, fasta_combined):
    """
    @summary: separate, FROGS artiificila combined OTU (replacing N by "-") and full length OTU.
    @param input_fasta : [str] Path to input fasta file
    @param fasta_full_length : [str] Path to output fasta file of full length OTU
    @param fasta_combined : [str] Path to output fasta file of artiificial combined OTU
    """

    FH_input = FastaIO(input_fasta)
    FH_FL = FastaIO(fasta_full_length, "wt")
    FH_AC = FastaIO(fasta_combined, "wt")

    nb = 0
    nb_combined = 0
    for record in FH_input:
        nb += 1
        if "N" in record.string:
            nb_combined +=1
            record.string.replace("N","-")
            FH_AC.write(record)
        else:
            FH_FL.write(record)
    FH_input.close()
    FH_FL.close()
    FH_AC.close()

    return nb,nb_combined

def split_fasta(fasta_file, tmp_files_manager, nb, out_list, log_file, by_nb_reads=False):
    """
    @summary: split fasta in nb_file and returne outfile list
    @param fasta_file: [str] Path to the fasta file to process.
    @param tmp_files_manager: [TmpFiles] The temporary file manager.
    @param nb : [int] The number of file to generate or number of reads per file.
    @param out_list : [list] List of output fasta file
    @param log_file : [srt] path to logfile
    @param by_nb_reads : [bool] choose to split input file in a fixed number of files or with a fixed number of reads per file
    """
    out_files = list()
    record_iter = FastaIO(fasta_file)
    for idx, record in enumerate(record_iter):
        out_file_idx = int(idx / nb ) if by_nb_reads else idx % nb
        if len(out_files) == 0 or not out_file_idx < len(out_files):
            if fasta_file.startswith(tmp_files_manager.prefix):
                new_out_file = tmp_files_manager.add( os.path.basename(fasta_file) + "_" + str(out_file_idx) , prefix="" )
            else: 
                new_out_file = tmp_files_manager.add( os.path.basename(fasta_file) + "_" + str(out_file_idx) )
            out_files.append({ 'file_path': new_out_file,
                               'file_handle': FastaIO(new_out_file, "wt"),
                               'nb_seq': 0
            })
            out_list.append( new_out_file )
        out_files[out_file_idx]['nb_seq'] += 1
        out_files[out_file_idx]['file_handle'].write( record )
    for out_file in out_files:
        out_file['file_handle'].close()

    # Log
    FH_log = Logger(log_file)
    FH_log.write("########################################################################################################\n")
    FH_log.write("# split " + fasta_file + " in smaller fasta files\n")
    FH_log.write("Results\n")
    if by_nb_reads:
        if out_files[-1]['nb_seq'] != nb:
            FH_log.write("\tGenerate " + str(len(out_list)-1) + " fasta files of " + str(nb) + " reads\n")
            FH_log.write("\tGenerate 1 fasta file of " + str(out_files[-1]['nb_seq']) + " reads\n")
        else:
            FH_log.write("\tGenerate " + str(len(out_list)) + " fasta files of " + str(nb) + " reads\n")
    else:
        for out_file in out_files:
            FH_log.write( "\tWrote " + str(out_file['nb_seq']) + " records to " + out_file['file_path'] + "\n" )
    FH_log.close()

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
    FH_summary_out = open( summary_file, "wt" )
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
        if cluster["metadata"]["blast_taxonomy"] is not None and len(cluster["metadata"]["blast_taxonomy"]) > 0:
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
            if sample_name not in samples_results:
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
                if cluster["metadata"]["blast_taxonomy"] is not None and len(cluster["metadata"]["blast_taxonomy"]) > 0:
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

def rdp_parallel_submission( function, inputs, outputs, logs, cpu_used, reference, memory):
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
            raise Exception("\n\n#ERROR : Error in sub-process execution.\n\n")

def process_multiple_needleall(reference, inputs_fasta, temp_sams, temp_logs, outputs, log_files, tempFiles_manager, debug):    
    for idx in range(len(inputs_fasta)):
        process_needleall(reference, inputs_fasta[idx], temp_sams[idx], temp_logs[idx], outputs[idx], log_files[idx], tempFiles_manager, debug)
     
def process_needleall(reference, input_fasta, temp_sam, temp_log, output, log_file, tmpFiles_manager, debug):
    """
    @summary: Launches NeedleAll on best blast refence.
    @param reference: [str] Path to the reference fasta file.
    @param input_fasta: [str] Path to the query fasta file to submit to NeedleAll.
    @param temp_sam: [str] Path to the NeedleAll sam results.
    @param temp_log: [str] Path to the NeedleAll log.
    @param output: [str] Path to blast-like tsv converted NeedleAll results.
    @param log_file: [str] Path to log.
    """
    
    # global alignment
    needleall_cmd = Needleall(reference, input_fasta, temp_sam, temp_log)
    needleall_cmd.submit(log_file)
    # convert needleall sam to tsv (blast like) 
    convert_to_tsv_cmd = NeedleallSam_to_tsv(temp_sam, reference, output)
    convert_to_tsv_cmd.submit(log_file)
    
    if not debug:
        tmpFiles_manager.delete(temp_sam)

def needleall_parallel_submission( function, reference, inputs_fasta, temp_sams, temp_logs, outputs, log_files, tmpFiles_manager, debug, cpu_used):

    processes = [{'process':None, 'inputs_fasta':list(), 'temp_sams':list(), 'temp_logs' :list(), 'outputs':list(), 'log_files':list()} for idx in range(cpu_used)]

    # Launch processes
    for idx in range(len(inputs_fasta)):
        process_idx = idx % cpu_used
        processes[process_idx]['inputs_fasta'].append(inputs_fasta[idx])
        processes[process_idx]['temp_sams'].append(temp_sams[idx])
        processes[process_idx]['temp_logs'].append(temp_logs[idx])
        processes[process_idx]['outputs'].append(outputs[idx])
        processes[process_idx]['log_files'].append(log_files[idx])

    for current_process in processes:
        if idx == 0:  # First process is threaded with parent job
            current_process['process'] = threading.Thread(target=function,
                                                          args=(reference, current_process['inputs_fasta'], current_process['temp_sams'], current_process['temp_logs'], current_process['outputs'], current_process['log_files'], tmpFiles_manager, debug))
        else:  # Others processes are processed on diffrerent CPU
            current_process['process'] = multiprocessing.Process(target=function,
                                                                 args=(reference, current_process['inputs_fasta'], current_process['temp_sams'], current_process['temp_logs'], current_process['outputs'], current_process['log_files'], tmpFiles_manager, debug))
        current_process['process'].start()
    # Wait processes end
    for current_process in processes:
        current_process['process'].join()
    # Check processes status
    for current_process in processes:
        if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
            raise Exception("\n\n#ERROR : Error in sub-process execution.\n\n")

###################################################################################################################
###                                              MAIN                                                           ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Taxonomic affiliation of each OTU's seed by RDPtools and BLAST.")
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used. [Default: %(default)s]")
    parser.add_argument( '-m', '--java-mem', type=int, default=2, help="Java memory allocation in Go. [Default: %(default)s]")
    parser.add_argument( '-t', '--taxonomy-ranks', nargs='*', default=["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"], help='The ordered ranks levels present in the reference databank. [Default: %(default)s]' )
    parser.add_argument('--rdp', default=False,  action='store_true',  help="Use RDP classifier to affiliate OTU")
    parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-r', '--reference', required=True, help='Preformated reference file (format: blast-indexed FASTA).')
    group_input.add_argument('-b', '--input-biom', required=True, help='BIOM file (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help="FASTA file of OTU's seed (format: FASTA).")
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-biom', default='affiliation_abundance.biom', help='BIOM file with added affiliation annotations from blast/needleall and/or RDPtools. [Default: %(default)s]')
    group_output.add_argument('-s', '--summary', default='affiliation_OTU.html', help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument('-l', '--log-file', default=sys.stdout, help='The list of commands executed.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    # Temporary files
    tmpFiles = TmpFiles( os.path.split(args.output_biom)[0] )
    fasta_full_length = tmpFiles.add(os.path.basename(args.input_fasta + "_FROGS_full_length"))
    fasta_combined = tmpFiles.add(os.path.basename(args.input_fasta + "_FROGS_combined"))
    # rdp tmp
    fasta_rdp_list = []
    rdp_out_list = []
    log_rdp_list = []
    # needle on FROGS_combined
    fasta_needleall_list = []
    blast_combined_list = []
    sam_needleall_list = []
    log_needleall_list = []
    needleall_tsv_out_list = []
    log_process_needl_list = []
    # aln Blast on FROGS full length
    blast_out_list = []
    log_blast_list = []
    # merge aln
    aln_out = ""

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
        nb_seq, nb_combined = extract_FROGS_combined(args.input_fasta, fasta_full_length, fasta_combined)
        if nb_seq == 0 : 
            raise Exception("\n\n#ERROR : Your input fasta file is empty!\n\n")
        Logger.static_write(args.log_file, "Nb seq : " + str(nb_seq) + "\n")
        if nb_combined > 0 :
            Logger.static_write(args.log_file, "\t with nb seq artificially combined :" + str(nb_combined) +"\n")
        Logger.static_write(args.log_file,"\n")

        if args.nb_cpus == 1 or nb_seq < 10:
            # kmer method affiliation
            # RDP
            if args.rdp:
                rdp_out_list.append( tmpFiles.add(os.path.basename(args.input_fasta) + ".rdp") )
                process_rdp( args.input_fasta, rdp_out_list[0], args.log_file, args.reference, args.java_mem)
            # alignment method affiliation
            # global alignment
            if nb_combined > 0 :
                # blast R1 and R2
                split_combined_R1 = tmpFiles.add(os.path.basename(fasta_combined) + "_R1.fasta") 
                split_combined_R2 = tmpFiles.add(os.path.basename(fasta_combined) + "_R2.fasta") 
                SplitOnTag(fasta_combined , 100*"N" , split_combined_R1 , split_combined_R2 ).submit(args.log_file)
                blast_combined_R1 = tmpFiles.add(os.path.basename(split_combined_R1) + ".blast") 
                blast_combined_R2 = tmpFiles.add(os.path.basename(split_combined_R2) + ".blast") 
                Blast(args.reference, split_combined_R1, blast_combined_R1, args.nb_cpus).submit(args.log_file)
                Blast(args.reference, split_combined_R2, blast_combined_R2, args.nb_cpus).submit(args.log_file)
                # reducing ref based on blast alignment
                reduced_ref = tmpFiles.add("reduced_" + os.path.basename(args.reference))
                log_reducre_ref  = tmpFiles.add("reduced_" + os.path.basename(args.reference) + ".log")
                Reduce_Ref(args.reference,blast_combined_R1,blast_combined_R2,reduced_ref, log_reducre_ref).submit(args.log_file)
                # needleall alignment and conversion in tsv (blast like)
                sam_needleall_list.append( tmpFiles.add( os.path.basename(fasta_combined) + ".needleall.sam" ) )
                log_needleall_list.append( tmpFiles.add( os.path.basename(fasta_combined) + ".needleall.log" ) )
                needleall_tsv_out_list.append( tmpFiles.add( os.path.basename(fasta_combined) + ".needleall.blast_like" ) )
                process_needleall(reduced_ref, fasta_combined, sam_needleall_list[0], log_needleall_list[0], needleall_tsv_out_list[0], args.log_file, tmpFiles, args.debug)
                
            # local alignment  
            if nb_seq - nb_combined > 0 :               
                #BLAST
                blast_out_list.append( tmpFiles.add(os.path.basename(fasta_full_length) + ".blast") )
                Blast(args.reference, fasta_full_length, blast_out_list[0], args.nb_cpus).submit(args.log_file)
        # parallelisation
        else:
            # kmer method affiliation
            # RDP
            if args.rdp:
                split_fasta(args.input_fasta, tmpFiles, max(1, int(args.nb_cpus/3)), fasta_rdp_list, args.log_file)
                rdp_out_list = [tmpFiles.add(os.path.basename(current_fasta) + ".rdp") for current_fasta in fasta_rdp_list]
                log_rdp_list = [tmpFiles.add(os.path.basename(current_fasta) + "_rdp.log") for current_fasta in fasta_rdp_list]
                rdp_parallel_submission( process_rdp, fasta_rdp_list, rdp_out_list, log_rdp_list, len(fasta_rdp_list), args.reference, args.java_mem )
            
            # alignment method affiliation
            # global alignment
            if nb_combined > 0:
                # blast R1 and R2
                split_combined_R1 = tmpFiles.add(os.path.basename(fasta_combined) + "_R1.fasta") 
                split_combined_R2 = tmpFiles.add(os.path.basename(fasta_combined) + "_R2.fasta") 
                SplitOnTag(fasta_combined , 100*"N" , split_combined_R1 , split_combined_R2 ).submit(args.log_file)
                blast_combined_R1 = tmpFiles.add(os.path.basename(split_combined_R1) + ".blast") 
                blast_combined_R2 = tmpFiles.add(os.path.basename(split_combined_R2) + ".blast") 
                Blast(args.reference, split_combined_R1, blast_combined_R1, args.nb_cpus).submit(args.log_file)
                Blast(args.reference, split_combined_R2, blast_combined_R2, args.nb_cpus).submit(args.log_file)
                # reducing ref based on blast alignment
                reduced_ref = tmpFiles.add("reduced_" + os.path.basename(args.reference))
                log_reducre_ref  = tmpFiles.add("reduced_" + os.path.basename(args.reference) + ".log")
                Reduce_Ref(args.reference,blast_combined_R1,blast_combined_R2,reduced_ref, log_reducre_ref).submit(args.log_file)
                #split input fasta for parallelisation
                
                split_fasta(fasta_combined, tmpFiles, 10, fasta_needleall_list, args.log_file, by_nb_reads=True)
                # needleall alignment and conversion in tsv (blast like)
                sam_needleall_list = [tmpFiles.add(os.path.basename(current_fasta) + ".needleall.sam", prefix="") for current_fasta in fasta_needleall_list ]
                log_needleall_list = [tmpFiles.add(os.path.basename(current_fasta) + ".needleall.log", prefix="") for current_fasta in fasta_needleall_list ]
                needleall_tsv_out_list = [tmpFiles.add(os.path.basename(current_fasta) + ".needleall.blast_like", prefix="") for current_fasta in fasta_needleall_list ]
                log_process_needl_list = [tmpFiles.add(os.path.basename(current_fasta) + ".process_needle.log", prefix="") for current_fasta in fasta_needleall_list ]
                needleall_parallel_submission( process_multiple_needleall, reduced_ref, fasta_needleall_list, sam_needleall_list, log_needleall_list, needleall_tsv_out_list, log_process_needl_list,tmpFiles, args.debug, min( len(fasta_needleall_list), args.nb_cpus ))
            # BLAST
            if nb_seq - nb_combined > 0 :
                blast_out_list.append(tmpFiles.add(os.path.basename(fasta_full_length) + ".blast") )
                log_blast_list.append( tmpFiles.add(os.path.basename(fasta_full_length) + "_blast.log") ) 
                Blast(args.reference, fasta_full_length, blast_out_list[0], args.nb_cpus).submit(log_blast_list[0])

            # Logs
            log_append_files(args.log_file, log_rdp_list + log_process_needl_list + log_blast_list)

        # Convert to output file
        AddAffiliation2Biom( args.reference, blast_out_list + needleall_tsv_out_list, rdp_out_list, args.input_biom, args.output_biom ).submit( args.log_file )
        summarise_results( args.summary, args.output_biom, args.taxonomy_ranks )

    finally:
        if not args.debug:
            tmpFiles.deleteAll()
