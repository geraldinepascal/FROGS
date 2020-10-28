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

__author__ = 'Maria Bernard - Sigenae team Jouy en Josas'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os,sys
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
def get_nb_seq( reads_file ):
    """
    @summary: Returns the number of sequences
    @param reads_file: [str] Path to the fasta/q file processed.
    @return: [int] The number of sequences.
    """
    FH_input = None
    if not is_gzip(reads_file):
        FH_input = open( reads_file, 'rt')
    else:
        FH_input = gzip.open( reads_file,'rt' )
    nb_line = 0
    for line in FH_input:
        nb_line += 1
    FH_input.close()

    format = "fastq" if FastqIO.is_valid(reads_file) else "fasta"
    nb_seq = nb_line/4 if format == "fastq" else nb_line/2
    return nb_seq

def checkQualityEncode ( input ):
    """
    @summary : check the quality Phred scale. Inspire on https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py
    @param input : [str] Path of fastq file
    @return : either "Q33" or "Q64"
    """

    RANGES = {
        'Q33': (33, 74),
        'Q64': (64, 104)
    }

    FH_in = FastqIO(input)
    # print(input)
    encoding = ""
    gmin, gmax  = 99, 0

    for record in FH_in : 
        qual = [ ord(c) for c in record.quality ]
        qmin = min(qual)
        qmax = max(qual)
        if qmin < gmin or qmax > gmax:
            gmin, gmax = min(qmin, gmin), max(qmax, gmax)
            valid_encodings = []
            for encoding, (emin, emax) in list(RANGES.items()):
                if gmin >= emin and gmax <= emax:
                    valid_encodings.append(encoding)
            if len(valid_encodings) == 1 and encoding == "" :
                encoding = valid_encodings[0]
                # break
    return encoding

def reverse_complement(sequence):
    """
    @summary : reverse complement DNA sequence
    @param sequence : [str] the sequence to reverse complement 
    @return the reversed complemented sequence
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'} 
    letters = list(sequence)
    letters = [basecomplement[base] for base in letters[::-1]]
    return ''.join(letters)

def splitSeq (input, format, tag, revcomp, out1, out2):
    """
    @summary : split input sequences on tag and ouput out1 and out2 or out1 interleaved pair sequences. Sequences ID will ends with _FROGS_split_part1 and _FROGS_split_part2
    @param input : [str] Path to fasta/q input sequences to split-output2
    @param format : [str] either fasta or fastq
    @param tag : [str] sequence tag to find in sequence to split them into 2 parts
    @param revcomp : [bool] To reverse complement R2 before combine
    @param out1 : [str] Path to part 1 splitted sequence or to interleaved parts of splitted sequences
    @param out2 : [str] Path to part 2 splitted sequence (optionnal)
    """

    FH_in = FastqIO(input) if format == "fastq" else FastaIO(input)
    FH_out1 = FastqIO(out1,"wt") if format == "fastq" else FastaIO(out1,"wt")
    if out2:
        FH_out2 = FastqIO(out2,"wt") if format == "fastq" else FastaIO(out2,"wt")

    # seq_id=[]

    for record in FH_in:
        record.id=record.id.replace("_FROGS_combined","")

        # checking for duplicate like that take too much time!!
        # if record.id in seq_id:
        #     raise Exception(record.id+" present multiple time in your input file")
        # else:
        #     seq_id.append(record.id)

        split_seq = record.string.split(tag)
        if len(split_seq) != 2:
            raise Exception("\n\n#ERROR : " + record.id + " of " + input + " can not be split into 2 pieces with the tag : "+ tag +"\n\n")

        if not record.description is None:
            R1_desc = record.description.split(";")[0].replace('R1_desc:','')
            R2_desc = record.description.split(";")[1].replace('R2_desc:','')
        else:
            R1_desc = None
            R2_desc = None

        if format == "fastq" :
            part1=Sequence(record.id+"_FROGS_split_part1", split_seq[0], R1_desc, record.quality[ : len(split_seq[0]) ] )
            if revcomp : 
                part2=Sequence(record.id+"_FROGS_split_part2", reverse_complement(split_seq[1]),R2_desc, record.quality[ :len(split_seq[0])+len(tag)-1:-1 ] )
            else :
                part2=Sequence(record.id+"_FROGS_split_part2", split_seq[1],R2_desc, record.quality[ len(split_seq[0])+len(tag) : ] )
        else:
            part1=Sequence(record.id+"_FROGS_split_part1", split_seq[0],R1_desc, None )
            if revcomp :
                part2=Sequence(record.id+"_FROGS_split_part2", reverse_complement(split_seq[1]),R2_desc, None )
            else :
                part2=Sequence(record.id+"_FROGS_split_part2", split_seq[1],R2_desc, None )

        FH_out1.write(part1)
        if out2:
            FH_out2.write(part2)
        else:
            FH_out1.write(part2)


def combineSeq(input1, input2, format, tag, revcomp, out):
    """
    @summary : combine 2 sequence (either in interleaved fasta/q file or in 2 files) with a specific tag (corresponding quality are automatically added in function of Phred scale encoding). Sequence ID will ends with "_FROGS_combined"
    @param input1  : [str] Path to sequence 5' or interleaved 5' / 3' sequences (fasta/q format)
    @param input2  : [str] Path to sequence 3' (fasta/q format) or None
    @param format  : [str] either "fasta" or "fastq"
    @param tag     : [str] the sequence tag to add between sequences
    @param revcomp : [bool] To reverse complement R2 before combine
    @param out     : [str] Path to fasta/q combined sequence output file
    """
    
    if format == "fastq":
        encode = checkQualityEncode( args.reads1 )
        badQualCode="!" if encode == "Q33" or encode == "" else "C"

    FH_in1 = FastqIO(input1) if format == "fastq" else FastaIO(input1)
    FH_in2 = None
    if input2 : 
        FH_in2 = FastqIO(input2) if format == "fastq" else FastaIO(input2)
    FH_out = FastqIO(out,"wt") if format == "fastq" else FastaIO(out,"wt")

    iter1 = FH_in1.__iter__()
    for record1 in iter1:
        record1.id = record1.id.replace("_FROGS_split_part1","")
        record1.desc = record1.description
        if record1.id.endswith(".1") or record1.id.endswith("/1"):
            record1.id=record1.id[:-2]

        if input2:
            record2 = FH_in2.next_seq()
        else : 
            record2 = next(iter1)

        record2.id = record2.id.replace("_FROGS_split_part2","")
        record2.desc = record2.description
        if record2.id.endswith(".2") or record2.id.endswith("/2"):
            record2.id=record2.id[:-2]

        if record1.id != record2.id:
            raise Exception ("\n\n#ERROR : Input files are not in correct order, starting with "+record1.id+" and "+record2.id+"\n\n")
        description = None
        if record1.desc != None and record2.desc != None :
            description = "R1_desc:"+record1.desc+";R2_desc="+record2.desc
        if format == "fastq" : 
            if revcomp:
                combined = Sequence(record1.id+"_FROGS_combined", \
                    record1.string+tag+reverse_complement(record2.string), \
                    description,\
                    record1.quality+badQualCode*len(tag)+record2.quality[::-1])
            else : 
                combined = Sequence(record1.id+"_FROGS_combined", \
                    record1.string+tag+record2.string, \
                    description,\
                    record1.quality+badQualCode*len(tag)+record2.quality)
        else : 
            combined = Sequence(record1.id+"_FROGS_combined", \
                record1.string+tag+record2.string, \
                description,\
                None)
        FH_out.write(combined)


def process( args ) :
    """
    @summary : Combine R1 and R2 reads in one sequence with a combine tag between them
               or 
               Split sequence on a split tag into 2 sequences
               or
               Split on split tag and combine on combine tag, so replacing one tag by an other in combined sequences
    """
    if args.combine_tag : 
        tmp_files = TmpFiles( os.path.split(args.combined_output)[0] )
    else :
        tmp_files = TmpFiles( os.path.split(args.split_output1)[0] )

    format = "fastq" if FastqIO.is_valid(args.reads1) else "fasta"
    if args.reads2:
        format2 = "fastq" if FastqIO.is_valid(args.reads1) else "fasta"
        if format != format2:
            raise Exception("\n\n#ERROR : Your reads1 and reads2 are not in the same format\n\n")

    try :
        # replace split tag by combine tag
        if args.split_tag and args.combine_tag :
            out_split1 = tmp_files.add( os.path.basename(args.reads1) + '_split1.'+format )
            out_split2 = tmp_files.add( os.path.basename(args.reads1) + '_split2.'+format )
            splitSeq(args.reads1, format, args.split_tag, False, out_split1, out_split2)
            combineSeq(out_split1, out_split2, format, args.combine_tag, False, args.combined_output)
        # split sequence on split tag
        elif args.split_tag :
            splitSeq(args.reads1, format, args.split_tag, True, args.split_output1, args.split_output2)
        # combine sequence with combine tag
        else :
            combineSeq(args.reads1, args.reads2, format, args.combine_tag, True, args.combined_output)
    finally:
        if not args.debug:
            tmp_files.deleteAll()


    
    
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( 
               description='Combine or split fasta/q sequences with/on tag sequence. If use both combine (C) and split (S) tag, the sequence will be split on S and combine with C (ie replacing S by C)')
    parser.add_argument( '-c', '--combine-tag', type=str, default=None, help='The sequence tag used to combine R1 and R2 [Default: %(default)s]' )
    parser.add_argument( '-s', '--split-tag', type=str, default=None, help='The sequence to find to split sequence in R1 and R2 [Default: %(default)s]' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '--debug', default=False, action='store_true', help="Keep temporary files to debug program." )
    group_input = parser.add_argument_group( 'Inputs' ) # Inputs
    group_input.add_argument( '--reads1', required=True, help='Path to read 1 or interleaved read1/read2 to combine\nOR\nPath to sequences to split' )
    group_input.add_argument( '--reads2', default=None, required=False, help='Path to read 2 file to combine (if not interleaved)' )
    group_output = parser.add_argument_group( 'Outputs' ) # Outputs
    group_output.add_argument( '--combined-output', help='Fasta/q output file of R1 and R2 reads combined with sequence tag <-c>. ')
    group_output.add_argument( '--split-output1', help='Fasta/q output file of :\n\t1) first part of splitted sequence on <-s> tag if --split-output2 specified.\n\t2) interleaved splitted sequences if --split-output2 not specified.')
    group_output.add_argument( '--split-output2', default=None, help='Fasta/q output file of second part of splitted sequence on <-s> tag. ')
    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')
    args = parser.parse_args()

    Logger.static_write(args.log_file, "## Application\nSoftware: " + os.path.basename(sys.argv[0]) + " (version: " + str(__version__) + ")\nCommand: " + " ".join(sys.argv) + "\n\n")

    ## check param
    #check tags, at least on of them, and at list one output
    if args.combine_tag is None and args.split_tag is None:
        raise Exception ("\n\n#ERROR : You need to provide at least one combine or split tag\n\n")
    if (not args.combined_output is None and not args.split_output1 is None) or (args.combined_output is None and args.split_output1 is None ):
        raise Exception ("\n\n#ERROR : You need to choose between combined output or split output\n\n")
    

    # check combine io
    if not args.reads2 is None and args.combine_tag is None:
        raise Exception("\n\n#ERROR : You provide 2 input files but no combine tag to combine them\n\n")

    if args.combine_tag is None and not args.combined_output is None :
        raise Exception ("\n\n#ERROR : You provide combined_output file name but no combined tag!\n\n")

    if not args.combine_tag is None and args.combined_output is None :
        raise Exception ("\n\n#ERROR : You provide combined_tag but no combined_output file name!\n\n")

    # check splitting io
    if args.split_output1 is None and not args.split_output2 is None:
        raise Exception("\n\n#ERROR : You can not provide split-output2 file name without split-output1 file name\n\n")

    process(args)  

