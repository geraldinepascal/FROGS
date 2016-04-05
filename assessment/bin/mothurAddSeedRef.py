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
__status__ = 'beta'


import re
import argparse
from sequenceIO import *

##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Add the original reference ID to each Mothur seed.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-i', '--input', required=True, help='Sequences file from mothur get.oturep (format: fasta).' )
    parser.add_argument( '-t', '--trimmed-reads', required=True, nargs="+", help='Reads after all sequence modifications like aln (format: fasta or fastq). It is used to find the ID of clusters centroids by exact comparison between OTU sequences and reads sequences.' )
    parser.add_argument( '-r', '--reads', required=True, nargs="+", help='Simulated reads used as input in mothur pipeline (format: fasta or fastq). These reads are used to retrieve simulation reference of the centroids. The link between centroids and reads is the sequence ID. The description of reads must contain the tag "reference=<REF_ID>".' )
    parser.add_argument( '-o', '--output', required=True, help='Output file (format: fasta).' )
    args = parser.parse_args()
    
    # Get observation sequences
    nb_observations = 0
    observation_ids_by_seq = dict()
    FH_seeds = FastaIO(args.input)
    for record in FH_seeds:
        nb_observations += 1
        if record.string not in observation_ids_by_seq:
            observation_ids_by_seq[record.string] = list()
        observation_ids_by_seq[record.string].append(record.id)
    FH_seeds.close()

    # Get centroids (the real centroid and indentical sequences) ID by observation
    observation_ids_by_centroid_id = dict()
    for file in args.trimmed_reads:
        FH_reads = SequenceFileReader.factory(file)
        for record in FH_reads:
            record_seq = record.string.replace("-", "").replace(".", "")
            if record_seq in observation_ids_by_seq:
                observation_ids_by_centroid_id[record.id] = observation_ids_by_seq[record_seq]
        FH_reads.close()

    # Get reference by observation
    reference_by_observation_id = dict()
    for file in args.reads:
        FH_reads = SequenceFileReader.factory(file)
        for record in FH_reads:
            if record.id in observation_ids_by_centroid_id:
                observation_ids = observation_ids_by_centroid_id[record.id]
                reference_id = re.search("reference=([^\s]+)", record.description).group(1)
                for current_obs_id in observation_ids:
                    if current_obs_id not in reference_by_observation_id:
                        reference_by_observation_id[current_obs_id] = reference_id
                    elif len(reference_by_observation_id[current_obs_id].split(",")) > len(reference_id.split(",")):
                        reference_by_observation_id[current_obs_id] = reference_id
        FH_reads.close()    
    if nb_observations != len(reference_by_observation_id):
        raise Exception("All the centroids sequences cannot be retrieved in reads files.")

    # Write seeds fasta with reference information
    FH_seeds = FastaIO(args.input)
    FH_seeds_with_ref = FastaIO(args.output, "w")
    for record in FH_seeds:
        record.description = "reference=" + reference_by_observation_id[record.id]
        FH_seeds_with_ref.write(record)
    FH_seeds.close()
    FH_seeds_with_ref.close()
