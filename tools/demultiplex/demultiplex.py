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

__author__ = 'Plateforme bioinformatique Toulouse and SIGENAE'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.2.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import os
import sys
import gzip
import time
import tarfile
import argparse

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


##################################################################################################################################################
#
# COMMAND LINES
#
##################################################################################################################################################
class Demultiplex(Cmd):
    """
    @summary : Demultiplex samples.
    """
    def __init__(self, R1_input_file, R2_input_file, barcode_file, mismatches, end, global_tmp_files, R1_output_files, R2_output_files, demultiplex_err_files1, demultiplex_err_files2, demultiplex_log):
        """
        @param R1_input_file : [str] Path to the R1 fastq file.
        @param R2_input_file : [str] Path to the R2 fastq file.
        @param barcode_file : [str] Path to barcodes and samples (one line by sample) description file. Line format : SAMPLE_NAME<TAB>BARCODE.
        @param mismatches : [int] Number of mismatches allowed
        @param end : [str] barcode ends ? forward : bol or reverse : eol (def bol)
        @param global_tmp_files : [str] Path for R1 and R2 files.
        @param R1_output_files : [list] Paths to the R1 fastq files (one by sample). User provides an empty list.
        @param R2_output_files : [list] Paths to the R2 fastq files (one by sample). User provides an empty list.
        @param demultiplex_err_files : [list] Paths to the files with ambiguous and unmatched reads. User provides an empty list.
        """
        
        tmp_files = TmpFiles( global_tmp_files.tmp_dir )
        
        tmp_folder = os.path.join( global_tmp_files.tmp_dir, global_tmp_files.prefix + "_tmp", tmp_files.prefix )
        global_tmp_files.dirs.append(tmp_folder)
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)
        self.samples_names = list()
        # Update output data
        FH_barcode = open( barcode_file )
        for line in FH_barcode:
            sample_name, barcode = line.strip().rsplit(None, 1)
            R1_output_files.append( os.path.join(tmp_folder, sample_name + '_R1.fastq') )
            global_tmp_files.files.append(os.path.join(tmp_folder, sample_name + '_R1.fastq') )
            if R2_input_file != None:
                R2_output_files.append( os.path.join(tmp_folder, sample_name + '_R2.fastq') )
                global_tmp_files.files.append(os.path.join(tmp_folder, sample_name + '_R2.fastq'))
            self.samples_names.append( sample_name.replace(' ', '_') )
        FH_barcode.close()
        self.R1_input_file = R1_input_file
        self.ambiguous = os.path.join(tmp_folder, 'ambiguous_R1.fastq')
        self.unmatched = os.path.join(tmp_folder, 'unmatched_R1.fastq')
        demultiplex_err_files1.extend( [self.ambiguous,self.unmatched] )
        global_tmp_files.files.extend( [self.ambiguous,self.unmatched] )
        if R2_input_file != None:
            demultiplex_err_files2.extend( [os.path.join(tmp_folder, 'ambiguous_R2.fastq'),os.path.join(tmp_folder, 'unmatched_R2.fastq') ])
            global_tmp_files.files.extend( [os.path.join(tmp_folder, 'ambiguous_R2.fastq'),os.path.join(tmp_folder, 'unmatched_R2.fastq') ])

        # Set class
        if R2_input_file != None:
            Cmd.__init__( self,
                          'splitbc.pl',
                          'Demultiplex reads.',
                          R1_input_file + ' ' + R2_input_file + ' --' + end + ' --bcfile ' + barcode_file + ' --mismatches ' + `mismatches` + ' --trim --no_adapt --prefix-r1 ' + os.path.join(tmp_folder, '%_R1.fastq') +\
                          ' --prefix-r2 ' + os.path.join(tmp_folder, '%_R2.fastq') + ' >> ' + demultiplex_log,
                          None )
        else:
            Cmd.__init__( self,
                          'splitbc.pl',
                          'Demultiplex reads.',
                          R1_input_file + ' --' + end + ' --bcfile ' + barcode_file + ' --mismatches ' + `mismatches` + ' --trim --no_adapt --prefix-r1 ' + os.path.join(tmp_folder, '%_R1.fastq') +\
                          ' >> ' + demultiplex_log,
                          None )
        
    def parser(self, log_file):
        """
        @summary : Parse the command results to add information in log_file.
        @log_file : [str] Path to the sample process log file.
        """
        # Parse output
        nb_seq_before = get_fastq_nb_seq(self.R1_input_file)
        nb_seq_unmatched = get_fastq_nb_seq(self.unmatched)
        nb_seq_ambiguous = get_fastq_nb_seq(self.ambiguous)
        # Write result
        FH_log = Logger( log_file )
        FH_log.write( 'Results :\n' )
        FH_log.write( '\tnb seq before demultiplexing : ' + str(nb_seq_before) + '\n' )
        FH_log.write( '\tnb seq after process matched : ' + str(nb_seq_before - nb_seq_unmatched) + '\n' )
        FH_log.write( '\tnb seq after process non-ambiguous : ' + str(nb_seq_before - nb_seq_unmatched - nb_seq_ambiguous) + '\n' )
        FH_log.close()

    def get_version(self):
        """
        @summary : Returns the program version number.
        @return : version number if this is possible, otherwise this method return 'unknown'.
        """
        return Cmd.get_version(self, 'stdout')


class Archive(Cmd):
    """
    @summary : Creates an archive with files.
    """
    def __init__(self, archived_files, archive_path):
        """
        @param archived_files: [list] Files added in final archive.
        @param archive_path: [str] Path to the new archive.
        """
        
        tmp_files=TmpFiles( os.path.dirname(archive_path) )
        tmp_folder = os.path.join( tmp_files.tmp_dir, tmp_files.prefix)
        tmp_files.dirs.append(tmp_folder)
        if not os.path.exists(tmp_folder):
            os.makedirs(tmp_folder)
            
        if len(archived_files) == 0:
            raise Exception( "At least one file must be add to the archive '" + archive_path + "'." )
    
        archived_basenames = list()
        for current in archived_files:
            if not os.path.dirname(current) == tmp_folder:
                os.rename(current, os.path.join(tmp_folder,os.path.basename(current)))
            tmp_files.files.append(os.path.join(tmp_folder,os.path.basename(current)))
            archived_basenames.append(os.path.basename(current))
                

        Cmd.__init__( self,
                      'tar',
                      'Archives files.',
                      '-zcf ' + archive_path + ' -C ' + tmp_folder + " " + " ".join(archived_basenames),
                      None )
        
        self.Files=tmp_files
        
        
    def parser(self,log_file):
        self.Files.deleteAll()
        

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def is_gzip( file ):
    """
    @return: [bool] True if the file is gziped.
    @param file : [str] Path to processed file.
    """
    is_gzip = None
    FH_input = gzip.open( file )
    try:
        FH_input.readline()
        is_gzip = True
    except:
        is_gzip = False
    finally:
        FH_input.close()
    return is_gzip

def split_barcode_file( barcode_file, barcodes_file_list, global_tmp_files ):
    """
    @summary: In case of double multiplexe, split barcode file in one forward and multiple reverse barcode files
    @param barcode_file: [str] Path to the input barcode file
    @param barcodes_file_list: [list] List of path to the ouput barcode files
    @param out_dir: [str] path to the output directory to write barcode files
    """
    out_dir = global_tmp_files.tmp_dir 
    barcode_input = open(barcode_file,"r")
    barcode_dict={}
    for l in barcode_input.readlines():
        [s,f,r]=l.strip().split()
        if not "forward_bc" in barcode_dict:
            barcode_dict["forward_bc"] = [f+"\t"+f]
        elif not f+"\t"+f in barcode_dict["forward_bc"]:
            barcode_dict["forward_bc"].append( f+"\t"+f)
        if not f+"_reverse_bc" in barcode_dict:
            barcode_dict[f+"_reverse_bc"] = [s+"\t"+r]
        else :
            barcode_dict[f+"_reverse_bc"].append(s+"\t"+r)

    f=barcode_dict.pop("forward_bc")
    barcodes_file_list.append(os.path.join(out_dir,"forward_bc"))
    global_tmp_files.files.append(os.path.join(out_dir,"forward_bc"))
    FH_out = open(os.path.join(out_dir,"forward_bc"),"w")
    FH_out.write("\n".join(f)+"\n")
    FH_out.close()

    for bc_file in barcode_dict:
        barcodes_file_list.append(os.path.join(out_dir,bc_file))
        global_tmp_files.files.append(os.path.join(out_dir,bc_file))
        FH_out = open(os.path.join(out_dir,bc_file),"w")
        FH_out.write("\n".join(barcode_dict[bc_file])+"\n")
        FH_out.close()

def get_fastq_nb_seq( fastq_file ):
    """
    @summary: Returns the number of sequences in fastq_file.
    @param fastq_file: [str] Path to the fastq file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(fastq_file):
        FH_input = open( fastq_file )
    else:
        FH_input = gzip.open( fastq_file )
    nb_line = 0
    for line in FH_input:
        nb_line += 1
    FH_input.close()
    nb_seq = nb_line/4
    return nb_seq

def concat_files(list_input, output_file):
    
    FH_out=open(output_file,"w")
    for f in list_input :
        FH_in = open(f)
        string=""
        i=0
        for line in FH_in:
            string+= line
            i+=1
            if i==2000 :
                FH_out.write(string)
                string=""
                i=0
        if i != 0:
            FH_out.write(string)
        FH_in.close()
    FH_out.close()

def summarise_results( summary_file, barcode_file, log_file ):
    """
    @summary: Writes one summary of results from several logs.
    @param summary_file: [str] The output file.
    @param log_files: [list] The list of path to log files (one log file by sample).
    """
    sample_dict=dict()
    FH_barcode= open(barcode_file)
    for line in FH_barcode:
        sample_dict[line.split()[0]]=0
    
    FH_summary = open(summary_file, "w")
    FH_summary.write( "#sample\tcount\n")
    FH_log = open(log_file,"r")
    sample_dict["unmatched"]=0
    sample_dict["ambiguous"]=0
    
    for line in FH_log.readlines():
        if line.startswith("Barcode") or  line.startswith("total") :
            pass
        else :
            l=line.replace('(','\t').split()
            if l[0] in sample_dict:
                sample_dict[l[0]] += int(l[1])
    
    for s in sample_dict:
        FH_summary.write(s + '\t' + str(sample_dict[s]) + '\n')
    FH_summary.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( 
        description='Split by samples the reads in function of inner barcode.'
    )
    parser.add_argument('-m', '--mismatches', type=int, default=0, help="Number of mismatches allowed in barcode")
    parser.add_argument('-e', '--end', type=str, default="bol", help="barcode is at the begining of the forward end (bol) or of the reverse (eol) or both (both). default forward (bol)")
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '--input-R1', required=True, help='The R1 sequence file with all samples (format: fastq).' )
    group_input.add_argument( '--input-R2', default=None, help='The R2 sequence file with all samples (format: fastq).' )
    group_input.add_argument( '--input-barcode', help='This file describes barcodes and samples (one line by sample). Line format : SAMPLE_NAME<TAB>BARCODE or SAMPLE_NAME<TAB>BARCODE_FW<TAB>BARCODE_RV.' )
    group_output = parser.add_argument_group( 'Outputs' )
    # Outputs
    group_output.add_argument( '--output-demultiplexed', default="demultiplexed_read.tar.gz", help='The tar file containing R1 files and R2 files for each sample (format: tar).' )
    group_output.add_argument( '--output-excluded', default="undemultiplexed_read.tar.gz", help='The tar file containing R1 files and R2 files not demultiplexed  (format: tar).' )
    group_output.add_argument( '-s', '--summary', default='summary.tsv', help='TSV file with summary of filters results  (format: TSV).')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()
    prevent_shell_injections(args)

    Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")

    # Process
    R1_files = list()
    R2_files = list()
    tmp_barcode_files = list()
    tmp_R1_files = list()
    tmp_R2_files = list()
    demultiplex_err_files1 = list()
    demultiplex_err_files2 = list()
    excluded_R1_file = os.path.join(os.path.split(args.output_demultiplexed)[0],os.path.basename(args.input_R1)+"_excluded_demult")
    if args.input_R2 != None :
        excluded_R2_file = os.path.join(os.path.split(args.output_demultiplexed)[0],os.path.basename(args.input_R2)+"_excluded_demult")
    uniq_id = str(time.time()) + "_" + str(os.getpid())
    
    tmp_files = TmpFiles( os.path.split(args.output_demultiplexed)[0] )
    demultiplex_log = tmp_files.add("Demult.log")
    tmp_folder=tmp_files.add_dir("tmp")
    os.mkdir(tmp_folder)
    
    sample_list=[]
    try:
        # Process
        if args.end == "bol" or args.end == "eol" :
            
            info="\n#Demultiplexing " + os.path.basename(args.input_R1)
            if args.input_R2 != None:
                info+= " and " + os.path.basename(args.input_R2) 
            info += " with " + os.path.basename(args.input_barcode) + " in " + args.end + " strand\n"
            Logger.static_write(args.log_file,info)
            Demultiplex(args.input_R1, args.input_R2, args.input_barcode, args.mismatches, args.end, tmp_files, R1_files, R2_files, demultiplex_err_files1,demultiplex_err_files2, demultiplex_log).submit( args.log_file )
        else:
            split_barcode_file(args.input_barcode, tmp_barcode_files, tmp_files)
            info="\n#Demultiplexing " + os.path.basename(args.input_R1)
            if args.input_R2 != None:
                info+= " and " + os.path.basename(args.input_R2) 
            info += " with " + os.path.basename(tmp_barcode_files[0]) + " in bol strand\n"
            Logger.static_write(args.log_file,info)
            Demultiplex(args.input_R1, args.input_R2, tmp_barcode_files[0], args.mismatches, "bol", tmp_files, tmp_R1_files, tmp_R2_files, demultiplex_err_files1,demultiplex_err_files2, demultiplex_log).submit( args.log_file )
            for idx,read1_file in enumerate(tmp_R1_files):
                bc = os.path.basename(read1_file).replace("_R1.fastq","")
                if os.path.join(tmp_files.tmp_dir,bc+"_reverse_bc") in tmp_barcode_files:
                    if os.stat(tmp_R1_files[idx]).st_size != 0 :
                        info="\n#Demultiplexing " + os.path.basename(tmp_R1_files[idx])
                        if args.input_R2 != None:
                            info+= " and " + os.path.basename(tmp_R2_files[idx])
                        info += " with " + bc+"_reverse_bc" + " in eol strand\n"
                        Logger.static_write(args.log_file,info)
                        if args.input_R2 != None:
                            Demultiplex(tmp_R1_files[idx], tmp_R2_files[idx], os.path.join(tmp_files.tmp_dir,bc+"_reverse_bc"), args.mismatches, "eol", tmp_files, R1_files, R2_files, demultiplex_err_files1, demultiplex_err_files2, demultiplex_log).submit( args.log_file )
                        else:
                            Demultiplex(tmp_R1_files[idx], None, os.path.join(tmp_files.tmp_dir,bc+"_reverse_bc"), args.mismatches, "eol", tmp_files, R1_files, R2_files, demultiplex_err_files1, demultiplex_err_files2, demultiplex_log).submit( args.log_file )
        
        Logger.static_write(args.log_file,"\n#Summarising result\n")
        summarise_results( args.summary, args.input_barcode, demultiplex_log )
        
        Logger.static_write(args.log_file,"\n#Concatenation of undemultiplexed files 1\n")
        concat_files(demultiplex_err_files1, excluded_R1_file )
        if len(R2_files) > 0:
            Logger.static_write(args.log_file,"\n#Concatenation of undemultiplexed files 2\n")
            concat_files(demultiplex_err_files2, excluded_R2_file )
            Logger.static_write(args.log_file,"\n#Archive demultiplexed R1 and R2 files\n")
            Archive(R1_files + R2_files, args.output_demultiplexed).submit( args.log_file )
            Logger.static_write(args.log_file,"\n#Archive undemultiplexed R1 and R2 files\n")
            Archive([excluded_R1_file,excluded_R2_file], args.output_excluded).submit( args.log_file )
        else:
            Logger.static_write(args.log_file,"\n#Archive demultiplexed files\n")
            Archive(R1_files, args.output_demultiplexed).submit( args.log_file )
            Logger.static_write(args.log_file,"\n#Archive undemultiplexed files\n")
            Archive([excluded_R1_file], args.output_excluded).submit( args.log_file )

    # Remove temporary files
    finally:
        if not args.debug:
            Logger.static_write(args.log_file,"\n#Removing temporary files\n")
            tmp_files.deleteAll()

