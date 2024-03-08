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

__author__ = 'Maria Bernard INRA - SIGENAE '
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '4.1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import sys
import json
import argparse

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


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class BIOM_sampling(Cmd):
    """
    @summary: Random sampling in each sample.
    """
    def __init__(self, in_biom, out_biom, nb_read, sampling_by_min, delete_samples):
        """
        @param in_biom: [str] Path to BIOM file.
        @param out_biom: [str] Path to output BIOM file.
        @param nb_read : [int] Number of reads per sample.
        @param sampling_by_min : [boolean] Sampling by the number of the smallest sample.
        """
        argument = ''
        if nb_read is not None:
            argument = " --nb-sampled " + str(nb_read)
        if sampling_by_min :
            argument = " --sampling-by-min "
        if delete_samples:
            argument += " --delete-samples"
        Cmd.__init__( self,
                      "biomTools.py",
                      "Random sampling in each sample.",
                      "sampling  --input-file " + in_biom + " --output-file " + out_biom + argument,
                      "--version" )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip()                      

class BIOM_FASTA_update(Cmd):
    """
    @summary: Converts BIOM file to TSV file.
    @note: taxonomyRDP seedID seedSequence blastSubject blastEvalue blastLength blastPercentCoverage blastPercentIdentity blastTaxonomy ASVname SommeCount sample_count
    """
    def __init__(self, in_biom, in_fasta, out_fasta, log):
        """
        @param in_biom: [str] Path to BIOM file.
        @param nb_read : [int] Number of reads per sample
        @param out_biom: [str] Path to output BIOM file.
        """
        Cmd.__init__( self,
                      'biomFastaUpdate.py',
                      'Update fasta file based on sequence in biom file',
                      "--input-biom " + in_biom + " --input-fasta " + in_fasta + " --output-file " + out_fasta + " --log " + log,
                      '--version' )

    def get_version(self):   
        return Cmd.get_version(self, 'stdout').strip() 

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_log(in_biom, num_reads, out_biom, log):
    sampling_by_min = False
    if num_reads is None:
        sampling_by_min = True
    FH_log=open(log,"wt")
    FH_log.write("#sample\tnb_asv_before\tnb_asv_after\n")
    initial_biom = BiomIO.from_json( in_biom )
    new_biom = BiomIO.from_json( out_biom )
    tot_seqs_before = 0
    tot_seqs_after = 0

    for sample_name in initial_biom.get_samples_names():
        if sample_name in new_biom.get_samples_names():
            nb_otu_before = len([ i for i in initial_biom.get_sample_obs(sample_name) if i >0 ])
            nb_seqs_before = sum([ i for i in initial_biom.get_sample_obs(sample_name) if i >0 ])
            tot_seqs_before += nb_seqs_before
            nb_otu_after = len([ i for i in new_biom.get_sample_obs(sample_name) if i > 0])
            tot_seqs_after += sum([ i for i in new_biom.get_sample_obs(sample_name) if i >0 ])
            if sampling_by_min is True or nb_seqs_before >= num_reads:
                FH_log.write("Sample name: "+sample_name+"\n\tnb initials ASV: "+str(nb_otu_before)+"\n\tnb normalised ASV: "+str(nb_otu_after)+"\n")
            else:
                FH_log.write("Below threshold sample: "+sample_name+"\n\tnb sequences: "+str(initial_biom.get_sample_count(sample_name))+"\n")
                FH_log.write("Sample name: "+sample_name+"\n\tnb initials ASV: "+str(nb_otu_before)+"\n\tnb normalised ASV: "+str(nb_otu_after)+"\n")
        else:
            tot_seqs_before += sum([ i for i in initial_biom.get_sample_obs(sample_name) if i >0 ])
            FH_log.write("Below threshold sample: "+sample_name+"\n\tnb sequences: "+str(initial_biom.get_sample_count(sample_name))+"\n")
            Logger.static_write(args.log_file,"WARNING: Deleted sample: "+str(sample_name) + " (Only " + str(initial_biom.get_sample_count(sample_name)) + " sequences).\n")

    nb_initial_otu=len(initial_biom.rows)
    nb_new_otu=len(new_biom.rows)
    nb_seqs_removed = tot_seqs_before - tot_seqs_after
    nb_otus_removed = nb_initial_otu - nb_new_otu
    FH_log.write("Sample name: all samples\n\tnb sequences kept: "+str(tot_seqs_after)+"\n\tnb sequences removed: "+str(nb_seqs_removed)+"\n\tnb ASV kept: "+str(nb_new_otu)+"\n\tnb ASV removed: "+str(nb_otus_removed)+"\n")
    FH_log.close()

def summarise_results( summary_file, is_delete_samples, num_reads, biom_subsample_log ):
    """
    @summary: Writes one summary of results from several logs.
    @param summary_file: [str] The output file.
    @param log_files: [list] The list of path to log files (one log file by sample).
    """
    # Get data
    # to summary ASVs number && abundances number              
    categories = ["Nb ASV before normalisation" ,"Nb ASV after normalisation" ]
    delete_categories = ['Nb sequences']
    series = list()
    deletes = list()
    get_sample_resuts( biom_subsample_log, series )
    i = 0
    n = len(series)
    while i < n :
        if series[i]["name"]=="all samples":
            summary = series.pop(i)["data"]
            n -= 1
        elif series[i]["name"].startswith('Below threshold sample:'):
            series[i]["name"] = series[i]["name"].replace('Below threshold sample:','')
            deletes.append(series.pop(i))
            n -= 1
        else:
            i += 1
    summary_info = {
       'abundance_kept' : summary[0],
       'abundance_removed' : summary[1],
       'nb_kept' : summary[2],
       'nb_removed' : summary[3]
    }
    # Write Deleted samples (nb sequences < ###NB_SEQS###)
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "normalisation_tpl.html") )
    FH_summary_out = open( summary_file, "wt" )
    if is_delete_samples:
        title = 'Deleted samples (nb sequences < ' + str(num_reads) + ')'
    else:
        title = 'Kept samples (nb sequences < ' + str(num_reads) + ')'

    for line in FH_summary_tpl:
        if "###DATA_CATEGORIES###" in line:
            line = line.replace( "###DATA_CATEGORIES###", json.dumps(categories) )
        elif "###DELETE_CATEGORIES###" in line:
            line = line.replace( "###DELETE_CATEGORIES###", json.dumps(delete_categories) )
        elif "###DATA_SERIES###" in line:
            line = line.replace( "###DATA_SERIES###", json.dumps(series) )
        elif "###DELETE_SERIES###" in line:
            line = line.replace( "###DELETE_SERIES###", json.dumps(deletes) )
        elif "###REMOVE_DATA###" in line:
            line = line.replace( "###REMOVE_DATA###", json.dumps(summary_info) )
        elif "###TITLE###" in line:
            line = line.replace( "###TITLE###", title )
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()

def get_sample_resuts( log_file, output_list ):
    """
    @summary: Returns the sample results (number of sequences after each filters).
    @param log_file: [str] Path to a log file.
    @return: [dict] The sample results. Format: {'name':SPL_NAME, 'data':[INI_NB_SEQ, NB_SEQ_AFTER_FILTER_1, NB_SEQ_AFTER_FILTER_2]}.
    """
    results = {
        'name': None,
        'data': list()
    }
    FH_input = open(log_file)
    for line in FH_input:
        if line.strip().startswith('Sample name: '):
            results['name'] = line.split(':')[1].strip()
        # case if informations about all samples
        elif line.strip().startswith('nb sequences kept:'):
            results['data'].append( int(line.split(':')[1].strip()) )
        elif line.strip().startswith('nb sequences removed:'):
            results['data'].append( int(line.split(':')[1].strip()) )
        elif line.strip().startswith('nb ASV kept'):
            results['data'].append( int(line.split(':')[1].strip()) )
        # case if informations about kept sample
        elif line.strip().startswith('nb initials ASV:'):
            results['data'].append( int(line.split(':')[1].strip()) )
        # case if informations about deleted sample
        elif line.strip().startswith('Below threshold sample:'):
            results['name'] = "Below threshold sample: " + line.split(':')[1].strip()
        # end of every cases and add all informations
        elif line.strip().startswith('nb sequences:') or\
         line.strip().startswith('nb normalised ASV:') or \
         line.strip().startswith('nb ASV removed:'):
            results['data'].append( int(line.split(':')[1].strip()) )
            output_list.append(results)
            results = {
                       'name': None,
                       'data': list()
                       }
    FH_input.close()

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Normalisation in BIOM by random sampling.")
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--debug', default=False, action='store_true', help="Keep temporary files to debug program. [Default: %(default)s]")
    parser.add_argument('--num-reads', type=int, help="Number of sampled sequences by sample.")
    parser.add_argument('--sampling-by-min', default=False, action='store_true', help='Sampling by the number of sequences of the smallest sample. [Default: %(default)s]' )
    parser.add_argument('--delete-samples', default=False, action='store_true', help='Delete samples that have a number of sequences below the selected filter. [Default: %(default)s]')
    
    
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-biom', required=True, help='Abundances file to normalise (format: BIOM).')
    group_input.add_argument('-f', '--input-fasta', required=True, help='Sequences file to normalise (format: FASTA).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', default='normalisation_abundance.biom', help='Normalised abundances (format: BIOM). [Default: %(default)s]')
    group_output.add_argument('--output-fasta', default='normalisation.fasta', help='Normalised sequences (format: FASTA). [Default: %(default)s]')
    group_output.add_argument('--summary-file', default='normalisation.html', help='The HTML file containing the graphs. [Default: %(default)s]')
    group_output.add_argument('--log-file', default=sys.stdout, help='The list of commands executed. [Default: stdout]')
    args = parser.parse_args()
    prevent_shell_injections(args)

    tmp_files = TmpFiles( os.path.split(args.output_biom)[0] )

    # Process
    try:
        Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")
        Logger.static_write(args.log_file,'Application start: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        
        Logger.static_write(args.log_file,'\n#Normalisation calculation\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        tmp_subsampling = tmp_files.add( 'tmp_biom_subsample.log' )

        if args.num_reads is None and not args.sampling_by_min:
            raise_exception( Exception('\n\n#ERROR : --sampling-by-min or --num-reads must be provided.\n\n'))

        BIOM_sampling(args.input_biom, args.output_biom, args.num_reads, args.sampling_by_min, args.delete_samples).submit(args.log_file)
        write_log(args.input_biom, args.num_reads, args.output_biom, tmp_subsampling)
        Logger.static_write(args.log_file,'\tend: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
        tmp_fastaUpdate = tmp_files.add( 'tmp_fasta_update.log' )
        BIOM_FASTA_update(args.output_biom, args.input_fasta, args.output_fasta, tmp_fastaUpdate).submit(args.log_file)
        Logger.static_write(args.log_file,'\n#Summarise\n\tstart: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
        summarise_results( args.summary_file, args.delete_samples, args.num_reads, tmp_subsampling )
        Logger.static_write(args.log_file,'\tend: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n\n' )
        Logger.static_write(args.log_file,'Application end: ' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + '\n' )
    # Remove temporary files
    finally:
        if not args.debug:
            tmp_files.deleteAll()
