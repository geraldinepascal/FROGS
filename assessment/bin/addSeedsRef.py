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
__version__ = '1.1.1'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'beta'

import re
import argparse
from frogsSequenceIO import FastaIO, SequenceFileReader
from frogsBiom import BiomIO


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Add the original reference ID to each UPARSE seed.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )
    group_input.add_argument( '-s', '--seeds-fasta', required=True, help='The sequence of the OTU seeds (format: fasta).' ) 
    group_input.add_argument( '-r', '--reads', required=True, nargs="+", help='The joined reads used in UPARSE pipeline (format: fasta or fastq). ID of the original sequence used to create read must be in reads description "reference=ID".' )
    # Outputs
    group_output = parser.add_argument_group( 'Outputs' )
    group_output.add_argument( '-a', '--annotated-fasta', required=True, help='The sequence of the OTU seeds with reference id in description (format: fasta).' )
    args = parser.parse_args()

    # Get observation sequences
    observation_id_by_seq = dict()
    FH_seeds = FastaIO(args.seeds_fasta)
    for record in FH_seeds:
        if record.string in observation_id_by_seq:
            raise Exception("The OTU '" + observation_id_by_seq[record.string] + "' and '" + record.id + "' have the same sequence.")
        observation_id_by_seq[record.string] = record.id.split(";size=")[0]
    FH_seeds.close()

    # Get centroids of observation
    reference_by_observation_id = dict()
    for file in args.reads:
        FH_reads = SequenceFileReader.factory(file)
        for record in FH_reads:
            if record.string in observation_id_by_seq:
                observation_id = observation_id_by_seq[record.string]
                reference_id = re.search("reference=([^\s]+)", record.description).group(1)
                if observation_id not in reference_by_observation_id:
                    reference_by_observation_id[observation_id] = reference_id
                elif len(reference_by_observation_id[observation_id].split(",")) > len(reference_id.split(",")):
                    reference_by_observation_id[observation_id] = reference_id
        FH_reads.close()
    if len(observation_id_by_seq) != len(reference_by_observation_id):
        missing = list()
        for seed_seq in observation_id_by_seq:
            if observation_id_by_seq[seed_seq] not in reference_by_observation_id:
                missing.append( observation_id_by_seq[seed_seq] )
        raise Exception("All the centroids sequences cannot be retrieved in reads files. Centroids without read: '" + "' '".join(missing) + "'.")

    # Write seeds fasta with reference information
    FH_seeds = FastaIO(args.seeds_fasta)
    FH_seeds_with_ref = FastaIO(args.annotated_fasta, "w")
    for record in FH_seeds:
        record.id = record.id.split(";size=")[0]
        record.description = "reference=" + reference_by_observation_id[record.id]
        FH_seeds_with_ref.write(record)
    FH_seeds.close()
    FH_seeds_with_ref.close()
