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
__email__ = 'support.genopole@toulouse.inra.fr'
__status__ = 'prod'

import argparse
from frogsSequenceIO import *


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Write groups of redundant sequences.' )
    parser.add_argument( '-i', '--input', type=str, required=True, help='Path to sequence file (format: fasta).' )
    parser.add_argument( '-o', '--output', type=str, required=True, help='Path to the output file (format: TSV). Each line is a redundant group (IDs of identical sequences).' )
    args = parser.parse_args()

    # Process
    #   Load IDs by sequence
    IDs_by_seq = dict()
    FH_databank = FastaIO(args.input)
    for record in FH_databank:
        if record.string not in IDs_by_seq:
            IDs_by_seq[record.string] = list()
        IDs_by_seq[record.string].append(record.id)
    FH_databank.close()
    #   Write redundant groups
    FH_output = open(args.output, "w")
    for sequence in IDs_by_seq:
        if len(IDs_by_seq[sequence]) > 1:
            FH_output.write( "\t".join(IDs_by_seq[sequence]) + "\n" )
    FH_output.close()
