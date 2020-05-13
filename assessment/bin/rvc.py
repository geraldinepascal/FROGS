#!/usr/bin/env python2.7
#
# Copyright (C) 2016 INRA
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

__author__ = 'Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'


import argparse
from frogsSequenceIO import *


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
complement_rules = {'A':'T','T':'A', 'G':'C','C':'G', 'U':'A', 'N':'N', 'Y':'R','R':'Y', 'S':'S', 'W':'W', 'K':'M','M':'K', 'B':'V','V':'B', 'D':'H','H':'D',
                    'a':'t','t':'a', 'g':'c','c':'g', 'u':'a', 'n':'n', 'y':'r','r':'y', 's':'s', 'w':'w', 'k':'m','m':'k', 'b':'v','v':'b', 'd':'h','h':'d'}

def rvc( in_file, out_file):
    """
    @summary: Writes the reversed complemented sequence.
    @param in_file: [str] Path to the input sequence file (fasta or fastq).
    @param out_file: [str] Path to the output file with reversed complemented sequences (the output format depend on input).
    """
    FH_in = SequenceFileReader.factory(in_file)
    FH_out = None
    if issubclass(FH_in.__class__, FastaIO):
        FH_out = FastaIO(out_file, "w")
    else:
        FH_out = FastqIO(out_file, "w")
    for record in FH_in:
        record.string = "".join([complement_rules[nt] for nt in record.string[::-1]])
        if record.quality is not None:
            record.quality = record.quality[::-1]
        FH_out.write(record)
    FH_in.close()
    FH_out.close()


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Reverse and complement sequences.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-i', '--input', required=True, help='Sequences file (format: fasta or fastq).' ) 
    parser.add_argument( '-o', '--output', required=True, help='Output file (format: same as input format).' )
    args = parser.parse_args()

    # Process
    rvc(args.input, args.output)
