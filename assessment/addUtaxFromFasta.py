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
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'

import re
import argparse
from sequenceIO import FastaIO
from biom import BiomIO


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Add taxonomy from UTAX result in BIOM file.")
    parser.add_argument( '-t', '--taxonomy-tag', default="taxonomy", help="The taxonomy tag in BIOM file. [Default: taxonomy]")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-f', '--input-fasta', required=True, help='Path to the sequence file outputed by UTAX (format: fasta).')
    group_input.add_argument('-b', '--input-biom', required=True, help='Path to the abundance file (format: BIOM).')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-biom', required=True, help='Path to the abundance file with taxonomy (format: BIOM).')
    args = parser.parse_args()

    # Process
    biom = BiomIO.from_json( args.input_biom )
    fasta = FastaIO( args.input_fasta )
    for record in fasta:
        # record.id example: Cluster_1;size=19714;tax=d:Bacteria(1.0000),p:"Proteobacteria"(0.9997),c:Alphaproteobacteria(0.9903),o:Rhodospirillales(0.9940),f:Acetobacteraceae(0.9887),g:Humitalea(0.9724);
        match = re.search("^([^\;]+)\;size\=\d+\;tax=(.+)$", record.id)
        if match is None:
            fasta.close()
            raise Exception("ID and taxonomy cannot be retrieved from '" + record.id + "'")
        record.id = match.group(1)
        record.description = match.group(2)
        biom.add_metadata( record.id, args.taxonomy_tag, record.description, "observation" )
    fasta.close()
    BiomIO.write( args.output_biom, biom )
