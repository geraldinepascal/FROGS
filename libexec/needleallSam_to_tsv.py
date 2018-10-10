#!/usr/bin/env python2.7
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

__author__ = 'Maria Bernard INRA - SIGENAE'
__copyright__ = 'Copyright (C) 2018 INRA'
__license__ = 'GNU General Public License'
__version__ = 'r3.0-v1.0'
__email__ = 'frogs@inra.fr'
__status__ = 'dev'

import os, sys
import argparse
import re

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
from frogsSequenceIO import *

###################################################################################################################
###                                                 FUNCTIONS                                                   ###
###################################################################################################################
def get_ref_length(in_db, out_dict):
	"""
	@summary: compute reference fasta length
	@param in_db : [str] path to reference fasta file
	@param out_dict : [dict] dictionnary to store reference sequence length
	"""
	FH_in = FastaIO(in_db)
	for record in FH_in:
		out_dict[record.id] = len(record.string.strip())

def process(params):
	"""
	@summary : convert Needleall Sam output in Blast like outfmt 6 format : 
				qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
	@param params : [Namespace] param.needle : needle input file, params.reference : needle refence fasta file and param.blast_like output file
	"""
	ref_length = dict()
	get_ref_length(params.reference,ref_length)
	OTU_list=list()
	temp_dict=dict()
	FH_in = open(params.needle)
	FH_out = open(params.blast_like,"w" )
	for line in FH_in:
		if line.startswith("@"):
			continue
		l = line.strip().split("\t")
		qseqid=l[0]
		sseqid=l[2]
		qstart=str(1)
		qend=str(len(l[9]))
		qlen=str(len(l[9]))

		bitscore=l[11].split(":")[-1]
		evalue=str(-1)

		cigar=l[5]
		cigarettes = re.findall('[\d]{0,}[A-Z]{1}', cigar)
		# external gap gives sstart and ssend
		if cigarettes[0].endswith("D"):
			sstart=str(int(cigarettes[0][:-1])+1)
			cigarettes.pop(0)
		else:
			sstart=str(1)
		if cigarettes[-1].endswith("D"):
			send=str(ref_length[sseqid]-int(cigarettes[-1][:-1]))
			cigarettes.pop(-1)
		else:
			send=str(ref_length[sseqid])
		# parse cigar, count gap, match and total length
		cig_dict = {"M":0, "D":0, "I":0}
		gapopen=0
		for cigarette in cigarettes:
			if cigarette.endswith("M"):
				cig_dict["M"]+= int(cigarette[:-1])
			elif cigarette.endswith("D"):
				cig_dict["D"]+= int(cigarette[:-1])
				gapopen += 1
			elif cigarette.endswith("I"):
				cig_dict["I"]+= int(cigarette[:-1])

		length = sum(cig_dict.values())

		# match / mismatch, and id%
		mismatch=l[12].split(":")[-1]

		# count N in query sequence, substract from length
		N = l[9].count("N")

		# SOLUTION 1 : nombre de match / longueur alignement
		# pident=str(round((cig_dict["M"]- int(mismatch))*100.00/length,2))

		# SOLUTION 2 : nombre de match / (longueur alignement - les 100 N)
		# pident=str(round((cig_dict["M"]- int(mismatch))*100.00/length,2))
		# length -= N

		# SOLUTION 3 : nombre de match / (longueur de la seed - les 100 N)
		pident=str(round((cig_dict["M"]- int(mismatch))*100.00/(len(l[9])-N),2))

		# SOLUTION 4 : nombre de match / (longueur de la seed )
		#pident=str(round((cig_dict["M"]- int(mismatch))*100.00/len(l[9]),2))

		# write sorted alignement by bitscore
		if not qseqid in OTU_list:
			if len(temp_dict) > 0:
				for a in sorted(temp_dict, key=lambda i: float(temp_dict[i]), reverse=True):
					FH_out.write(a)
				temp_dict.clear()
			OTU_list.append(qseqid)
		aln=qseqid+"\t"+sseqid+"\t"+pident+"\t"+str(length)+"\t"+mismatch+"\t"+str(gapopen)+"\t"+qstart+"\t"+qend+"\t"+sstart+"\t"+send+"\t"+evalue+"\t"+bitscore+"\t"+qlen+"\n"
		temp_dict[aln]=bitscore
	
	for a in sorted(temp_dict, key=lambda i: float(temp_dict[i]), reverse=True):
		FH_out.write(a)



###################################################################################################################
###                                              MAIN                                                           ###
###################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Convert Needle Sam output format to Blast outfmt 6 like format")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
	# Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-n', '--needle', required=True, help='Needle Sam output file')
    group_input.add_argument('-r', '--reference', required=True, help='Reference fasta file')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-b', '--blast-like', required=True, help='Blast like tsv format. qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen')

    args = parser.parse_args()
    prevent_shell_injections(args)

    process(args)