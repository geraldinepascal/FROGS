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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse AND Maria Bernard - SIGENAE Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.6.1'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import re
import sys
import time
import operator
import argparse
import threading
import multiprocessing

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
def getSamplesFromDescFile( samples_desc_file ):
    """
    @summary: Returns the list of samples names and the list of their sequences files.
    @param samples_desc_file: [str]
    @return: [list] The first element is the list of samples names. The second is the list of the sequences files.
    """
    sequence_files = list()
    samples_names = list()
    FH_samples_desc_file = open( samples_desc_file, 'rt' )
    for line in FH_samples_desc_file:
        if not line.startswith("#"):
            sample_path, sample_name = line.strip().split("\t", 1)
            sequence_files.append( sample_path )
            samples_names.append( sample_name )
    FH_samples_desc_file.close()
    return samples_names, sequence_files


def splitFileBySample(input_count, fasta, size_separator, working_dir, process_prefix):
    """
    @summary : split fasta file by sample with abundance indicated in id thanks to size_sep
    @param input_count : [str] Path to count tsv file
    @param fasta : [str] Path to fasta file
    @param size_separator : [str] The size separator in the sequences ID if the sequences are already dereplicated in each sample.
    @param working_dir: [str] The output directory.
    @param process_prefix: [str] The prefix for tmp files.
    """
    sequence_files = list()
    samples_names = list()

    # count file saving
    dict_count = dict()
    FH_count = open(input_count, 'rt')
    for line in FH_count:
        if line.startswith("#id"):
            samples_names = line.strip().split("\t")[1:]
            sequence_files = [ os.path.join(working_dir, process_prefix + "tmp_sample_" + s + '.fasta') for s in samples_names ]
        else:
            seq_id = line.split("\t")[0]
            dict_count[seq_id] = { samples_names[idx] : int(count) for idx,count in enumerate(line.strip().split()[1:]) }
    FH_count.close()

    # split fasta file
    print(fasta)
    FH_in = FastaIO(fasta)
    FH_out = dict()

    for record in FH_in:
        record.id = record.id.split(size_separator)[0]
        samples = [ s for s in dict_count[record.id] if dict_count[record.id][s] > 0 ]
        for s in samples:
            if not s in FH_out :
                FH_out[s] = FastaIO(os.path.join(working_dir, process_prefix + "tmp_sample_" + s + '.fasta'), 'at')
            sample_record = Sequence(record.id, record.string, record.description)
            sample_record.id = sample_record.id + size_separator + str(dict_count[sample_record.id][s])
            FH_out[s].write(sample_record)

    return samples_names, sequence_files

def splitMultiplesFilesByLength(files, samples_names, working_dir, process_prefix):
    """
    @summary: Split by length in all files in a list of files.
    @param files: [list] Files to split.
    @param sample_name: [list] Sample name for each file.
    @param working_dir: [str] The output directory.
    @param process_prefix: [str] The prefix for tmp files.
    """
    for idx, current_file in enumerate(files):
        splitFileByLength(current_file, samples_names[idx], working_dir, process_prefix)

def splitFileByLength(initial_file, sample_name, working_dir, process_prefix):
    """
    @summary: Write each sequence in file in one file by length.
    @param initial_file: [str] The file to split.
    @param sample_name: [str] The sample name for the file.
    @param working_dir: [str] The output directory.
    @param process_prefix: [str] The prefix for tmp files.
    """
    splits_FH = dict()

    reader = FastaIO(initial_file)
    for record in reader:
        if record.description is not None:
            record.description = record.description + " sample_name=" + sample_name
        else:
            record.description = " sample_name=" + sample_name
        seq_length = len(record.string)
        # if seq_length not in splits_FH:
            # splits_FH[seq_length] = FastaIO(os.path.join(working_dir, process_prefix + "tmp_" + os.path.basename(initial_file) + '_length_' + str(seq_length) + '.fasta'), 'at')

        if seq_length not in splits_FH:
            seq_file = os.path.join(working_dir, process_prefix + "tmp_" + os.path.basename(initial_file) + '_length_' + str(seq_length) + '.fasta')
            try:
                splits_FH[seq_length] = FastaIO(seq_file, 'at')
            except OSError:
                # remove randomly a file and close it
                _, file_obj = splits_FH.popitem() 
                file_obj.close()
                splits_FH[seq_length] = FastaIO(seq_file, 'at')

        splits_FH[seq_length].write( record )

    for length in list(splits_FH.keys()):
        splits_FH[length].close()

def dereplicateMultiplesLengths( working_dir, lengths, samples_names, process_prefix, size_separator ):
    """
    @summary : Dereplicates all sequences with length in a list of lengths.
    @param working_dir : [str] The directory with samples files splitted by sequence length.
    @param lengths : [list] Only sequences with these lengths are dereplicated.
    @param samples_names : [list] All samples names.
    @param process_prefix : [str] The prefix for tmp files.
    @param size_separator : [str] The size separator in the sequences ID if the sequences are already dereplicated in each sample.
    """
    for current_length in lengths:
        dereplicateLength( working_dir, current_length, samples_names, process_prefix, size_separator )

def dereplicateLength( working_dir, length, all_samples, process_prefix, size_separator=None ):
    """
    @summary: Dereplicates all sequences with length equals to "length". This method produce a fasta file and a count file.
    @param working_dir: [str] The directory with samples files splitted by sequence length.
    @param length: [str] Only sequences with this length are dereplicated.
    @param all_samples: [list] All samples names.
    @param process_prefix: [str] The prefix for tmp files.
    @param size_separator: [str] The size separator in the sequences ID if the sequences are already dereplicated in each sample.
    """
    records_list = list()
    derep_fasta = os.path.join(working_dir, process_prefix + "tmp_derep_" + str(length) + ".fasta")
    count_file = os.path.join(working_dir, process_prefix + "tmp_count_" + str(length) + ".tsv")

    # List records
    for item in sorted(os.listdir(working_dir )):
        if os.path.isfile(os.path.join(working_dir , item)) and re.match(process_prefix + 'tmp_.+_length_' + str(length) + '.fasta', item):
            reader = FastaIO(os.path.join(working_dir , item))
            for record in reader:
                records_list.append({"id":record.id,"desc":record.description,"seq":record.string})

    # Sort
    records_list.sort(key=operator.itemgetter('seq'))

    # Dereplicate
    FH_derep = FastaIO(derep_fasta, "wt")
    FH_count = open(count_file, "wt")
    prev_record = None
    prev_seq = None
    for idx, record in enumerate(records_list):
        record['desc'], sample_name = record['desc'].rsplit('sample_name=')
        if record['desc'].strip() == "":
            record['desc'] = None
        count = 1
        if size_separator is not None:
            record['id'], count = record['id'].rsplit(size_separator)
            if count.endswith(';'):
                count = count[:-1]
            count = int(count)
        if record['seq'] == prev_seq:
            if sample_name not in prev_record['count']:
                prev_record['count'][sample_name] = 0
            prev_record['count'][sample_name] += count
            records_list[idx] = None
        else:
            if prev_record is not None:
                FH_count.write( prev_record['id'] + "\t" + "\t".join([str(prev_record['count'][sample]) for sample in all_samples]) + "\n" )
                new_id = prev_record['id'] + ";size=" + str(sum(prev_record['count'].values()))
                FH_derep.write( Sequence(new_id, prev_record['seq'], prev_record['desc']) )
                prev_record['count'] = None
            prev_record = record
            prev_record['count'] = {name:0 for name in all_samples}
            prev_record['count'][sample_name] = count
            prev_seq = record['seq']
    if prev_record is not None:
        FH_count.write( prev_record['id'] + "\t" + "\t".join([str(prev_record['count'][sample]) for sample in all_samples]) + "\n" )
        new_id = prev_record['id'] + ";size=" + str(sum(prev_record['count'].values()))
        FH_derep.write( Sequence(new_id, prev_record['seq'], prev_record['desc']) )
    FH_derep.close()
    FH_count.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( 
               description='Remove duplicated sequences in multi-samples experimentation. Sequences using a strictly full-length matching.\n' +
                           'The number of sequences represented by a unique sequence is add in the end of ID.\n' +
                           'The count by sample is produce in "count-file".', 
               usage='''derepSamples.py [-h] [-v]
          [-s SIZE_SEPARATOR] [-p NB_CPUS]
          --samples-ref PATH | --sequences-files PATH [PATH ...]
                                [--samples-names NAMES [NAMES ...]]
          [-d DEREPLICATED_FILE] [-c COUNT_FILE]
''')
    parser.add_argument( '-s', '--size-separator', default=None, help='The size separator in the sequences ID if the "sequences-files" are already dereplicated in each sample.' )
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help='The maximum number of CPUs used. [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' ) 
    group_input.add_argument( '--samples-ref', help='Path to the file containing the link between samples names and sequences files (format: TSV). Each line in this file describes a sample: sequence_file_path<TAB>sample_name.' )
    group_input.add_argument( '--sequences-files', nargs='+', help='The sequence file for each sample (format: fasta). This parameter is mutally exclusive with the samples-ref parameter.' )
    group_input.add_argument( '--samples-names', nargs='+', help='The sample name for each sequences-files.' )
    group_input.add_argument( '--samples-count', default=None, help='The samples count TSV file if already available. This parameter is mutually exclusive with the samples-names parameter')
    # Outputs 
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-d', '--dereplicated-file', default='dereplication.fasta', help='Fasta file with unique sequences. Each sequence has an ID ended with the number of initial sequences represented (example : ">a0101;size=10"). [Default: %(default)s]')
    group_output.add_argument( '-c', '--count-file', default='count.tsv', help='TSV file with count by sample for each unique sequence (example with 3 samples : "a0101<TAB>5<TAB>8<TAB>0"). [Default: %(default)s]')
    args = parser.parse_args()

    if args.sequences_files is None and args.samples_ref is None:
        raise parser.error( "\n\n#ERROR : The parameter '--samples-ref' or '--sequences-files' is required.\n\n" )
    elif args.sequences_files is not None and args.samples_ref is not None:
        raise parser.error( "\n\n#ERROR : The parameters '--samples-ref' and '--sequences-files' are mutually exclusive.\n\n" )
    elif args.samples_names is not None and args.samples_count is not None:
        raise parser.error( "\n\n#ERROR : The parameters '--samples-names' and '--input-count' are mutually exclusive.\n\n" )
    elif args.sequences_files is not None and len(args.sequences_files) != 1 and args.samples_count:
        raise parser.error("\n\n#ERROR : The parameter '--input-count' must correspond to only one fasta file provided with '--sequences-files'\n\n")

    # Process
    working_dir = os.path.dirname(os.path.abspath(args.dereplicated_file))
    process_prefix = str(time.time()) + "_" + str(os.getpid()) + "_"
    lengths = dict()

    try:
        sequences_files = list()
        # Samples names
        if args.sequences_files is None:
            args.samples_names, sequences_files = getSamplesFromDescFile( args.samples_ref )
        elif args.samples_names is None and args.samples_count is None :
            args.samples_names = [os.path.splitext(os.path.basename(current_file))[0] for current_file in args.sequences_files]
            sequences_files = args.sequences_files
        elif args.samples_count:
            if args.size_separator is None:
                args.size_separator = ";size="
            args.samples_names, sequences_files = splitFileBySample(args.samples_count, args.sequences_files[0], args.size_separator, working_dir, process_prefix)

        # Split by length
        nb_processes_used = min( len(sequences_files), args.nb_cpus )
        if nb_processes_used == 1:
            splitMultiplesFilesByLength( sequences_files, args.samples_names, working_dir, process_prefix )
        else:
            processes = [{'process':None, 'files':[], 'samples_names':[]} for idx in range(nb_processes_used)]
            # Set processes
            for idx in range(len(sequences_files)):
                process_idx = idx % nb_processes_used
                processes[process_idx]['files'].append( sequences_files[idx] )
                processes[process_idx]['samples_names'].append( args.samples_names[idx] )
            # Launch processes
            for idx, current_process in enumerate(processes):
                if idx == 0: # First process is threaded with parent job
                    current_process['process'] = threading.Thread( target=splitMultiplesFilesByLength, 
                                                                   args=(current_process['files'], current_process['samples_names'], working_dir, process_prefix) )
                else: # Others processes are processed on diffrerent CPU
                    current_process['process'] = multiprocessing.Process( target=splitMultiplesFilesByLength, 
                                                                          args=(current_process['files'], current_process['samples_names'], working_dir, process_prefix) )
                current_process['process'].start()
            # Wait processes end
            for current_process in processes:
                current_process['process'].join()
            # Check processes status
            for current_process in processes:
                if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
                    raise Exception( "\n\n#ERROR : Error in sub-process execution.\n" )

        # List length
        for item in sorted(os.listdir(working_dir)):
            if os.path.isfile(os.path.join(working_dir, item)) and re.match(process_prefix + 'tmp_.+_length_\d+.fasta', item):
                current_length = item.rsplit('_length_', 1)[1][:-6]
                lengths[current_length] = True

        # Dereplicate by length
        nb_processes_used = min( len(list(lengths.keys())), args.nb_cpus )
        if nb_processes_used == 1:
            dereplicateMultiplesLengths( working_dir, list(lengths.keys()), args.samples_names, process_prefix, args.size_separator )
        else:
            processes = [{'process':None, 'length':[]} for idx in range(nb_processes_used)]
            # Set processes
            for idx, current_length in enumerate(lengths.keys()):
                process_idx = idx % nb_processes_used
                processes[process_idx]['length'].append(current_length)
            # Launch processes
            for current_process in processes:
                if idx == 0: # First process is threaded with parent job
                    current_process['process'] = threading.Thread( target=dereplicateMultiplesLengths, 
                                                                   args=(working_dir, current_process['length'], args.samples_names, process_prefix, args.size_separator) )
                else: # Others processes are processed on diffrerent CPU
                    current_process['process'] = multiprocessing.Process( target=dereplicateMultiplesLengths, 
                                                                   args=(working_dir, current_process['length'], args.samples_names, process_prefix, args.size_separator) )
                current_process['process'].start()
            # Wait processes end
            for current_process in processes:
                current_process['process'].join()
            # Check processes status
            for current_process in processes:
                if issubclass(current_process['process'].__class__, multiprocessing.Process) and current_process['process'].exitcode != 0:
                    raise Exception( "\n\n#ERROR : Error in sub-process execution.\n\n" )

        # Merge results
        FH_derepli = open(args.dereplicated_file, "wt")
        FH_count = open(args.count_file, "wt")
        FH_count.write( "#id\t" + "\t".join(args.samples_names) + "\n" )
        for item in sorted(os.listdir(working_dir)):
            if os.path.isfile(os.path.join(working_dir, item)):
                if re.match(process_prefix + 'tmp_count_\d+.tsv', item):
                    FH_count_part = open(os.path.join(working_dir, item), 'rt')
                    for line in FH_count_part:
                        FH_count.write( line )
                    FH_count_part.close()
                elif re.match(process_prefix + 'tmp_derep_\d+.fasta', item):
                    FH_derepli_part = open(os.path.join(working_dir, item), 'rt')
                    for line in FH_derepli_part:
                        FH_derepli.write( line )
                    FH_derepli_part.close()
        FH_count.close()
        FH_derepli.close()
    finally:
        # Delete tmp files
        if not args.debug:
            for item in sorted(os.listdir(working_dir)):
                if os.path.isfile(os.path.join(working_dir, item)) and (re.match(process_prefix + 'tmp_.+_length_\d+.fasta', item) or re.match(process_prefix + 'tmp_count_\d+.tsv', item) or re.match(process_prefix + 'tmp_sample_.+.fasta', item) or re.match(process_prefix + 'tmp_derep_\d+.fasta', item)):
                    os.remove(os.path.join(working_dir, item))
