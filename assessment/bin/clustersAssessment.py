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
import os
import argparse
from sequenceIO import *
from biom import Biom, BiomIO


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def samples_from_dir(directory_path, sample_name_sep):
    """
    @summary: Parse the directory and return by sample the path to the sequence file.
    @param directory_path: [str] Path to the samples directory.
    @param sample_name_sep: [str] Separator between sample name and the rest of the filename.
    @returns: [dict] By sample name the path of the sequence file.
    @warning: The directory must only contains the samples files.
    """
    samples = dict()
    for filename in os.listdir(directory_path):
        if os.path.isfile(os.path.join(directory_path, filename)):
            sample_name = filename.split(".")[0].split(sample_name_sep)[0]
            samples[sample_name] = os.path.join(directory_path, filename)
    return samples


def get_ref_after_simulation(samples):
    """
    @summary: Returns the IDs of initial sequences presents after simulation.
    @param samples: [dict] By sample name the path of the sequence file
    @returns: [list] The first element is a dict of reference present in dataset.
              The second element is a dict of reference present by sample.
    """
    references = dict()
    references_by_sample = dict()
    
    for sample_name in samples:
        expected_reference = dict()
        FH_sample = SequenceFileReader.factory(samples[sample_name])
        for record in FH_sample:
            origin = re.search("reference=([^\s]+)", record.description).group(1)
            if "," not in origin:
                expected_reference[origin] = 1
                references[origin] = 1
        references_by_sample[sample_name] = expected_reference
        FH_sample.close()
    return references, references_by_sample


def get_uniq(fasta, references_by_sample):
    """
    
    """
    nb_uniq = 0
    nb_uniq_by_sample = dict()

    # Get sequences
    sequence_by_id = dict()
    fh_sequences = FastaIO( fasta )
    for record in fh_sequences:
        sequence_by_id[record.id] = record.string
    fh_sequences.close()

    # Get nb uniq global
    nb_uniq = len(set([sequence_by_id[id] for id in sequence_by_id]))

    # Get nb uniq by sample
    for sample_name in references_by_sample:
        sequences = dict()
        for reference_id in references_by_sample[sample_name]:
            sequences[sequence_by_id[reference_id]] = 1
        nb_uniq_by_sample[sample_name] = len(sequences)

    return nb_uniq, nb_uniq_by_sample


def get_global_retrieved(reference_by_obs_id, references):
    expected_retrieved_already_processed = dict()
    retrieved_already_processed = dict()
    nb_detected = 0
    nb_splits = 0
    for obs_id in reference_by_obs_id:
        nb_detected += 1
        if not "," in reference_by_obs_id[obs_id]: # Is not a chimera
            ref_id = reference_by_obs_id[obs_id]
            if ref_id in retrieved_already_processed:
                nb_splits += 1
            else:
                retrieved_already_processed[ref_id] = 1
                if ref_id in references:
                    expected_retrieved_already_processed[ref_id] = 1
    
    return {
        "expected_retrieved": len(expected_retrieved_already_processed),
        "retrieved": len(retrieved_already_processed),
        "detected": nb_detected,
        "splits": nb_splits
    }


def get_retrieved_by_sample( biom_file, reference_by_obs_id ):
    counts_by_sample = dict()
    biom = BiomIO.from_json( biom_file )
    for sample_name in biom.get_samples_names():
        nb_detected = 0
        nb_expected_retrieved = 0
        retrieved_already_processed = dict()
        for obs in biom.get_observations_by_sample( sample_name ):
            nb_detected += 1
            if not "," in reference_by_obs_id[obs['id']]: # Is not a chimera
                ref_id = reference_by_obs_id[obs['id']]
                if ref_id not in retrieved_already_processed: # Not already processed this reference
                    retrieved_already_processed[ref_id] = 1
                    if ref_id in references_by_sample[sample_name]:
                        nb_expected_retrieved += 1
        counts_by_sample[sample_name] = {
            "detected": nb_detected,
            "retrieved": len(retrieved_already_processed),
            "expected_retrieved": nb_expected_retrieved
        }
    return counts_by_sample


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser( description='Count by sample # OTU, # real OTU and # real splitted OTU.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-f', '--fasta', required=True, help='Sequences file (format: fasta). ID of the original sequence used to create seed must be in seed description "reference=ID".' ) 
    parser.add_argument( '-b', '--biom', required=True, help='Abundance file (format: BIOM).' )
    parser.add_argument( '-o', '--origin', required=True, help='Sequence file provided to the simulation workflow (format: fasta).' )
    parser.add_argument( '-r', '--reads-dir', required=True, help='Path to the directory with one sequence file by sample (reads produced by simulation).' )
    parser.add_argument( '-s', '--sample-sep', default="_", help='Separator between sample name and the rest of the filename. [Default: _]' )
    args = parser.parse_args()
    
    # Get samples files
    samples = samples_from_dir(args.reads_dir, args.sample_sep)
    
    # Get expected reference
    references, references_by_sample = get_ref_after_simulation(samples)

    # Get uniq reference
    nb_uniq, nb_uniq_by_sample = get_uniq(args.origin, references_by_sample)

    # Get OTU in results
    reference_by_obs_id = dict()
    fh_sequences = FastaIO( args.fasta )
    for record in fh_sequences:
        if ";size=" in record.id:
            record.id = record.id.split(";size=", 1)[0]
        reference_by_obs_id[record.id] = re.search("reference=([^\s]+)", record.description).group(1)
    fh_sequences.close()
    counts = get_global_retrieved(reference_by_obs_id, references)
    counts_by_sample = get_retrieved_by_sample( args.biom, reference_by_obs_id )

    # Output
    print "#After_simu\tDictincts_after_simu\tExpected_retrieved\tRetrieved\tDetected\tSplits"
    print "\t".join([ str(len(references)),
                      str(nb_uniq),
                      str(counts["expected_retrieved"]),
                      str(counts["retrieved"]),
                      str(counts["detected"]),
                      str(counts["splits"]) ])
    print ""
    
    print "#Sample\tAfter_simu\tDistincts_after_simu\tExpected_retrieved\tRetrieved\tDetected"
    for sample in samples:
        print "\t".join([ sample, 
                          str(len(references_by_sample[sample])),
                          str(nb_uniq_by_sample[sample]),
                          str(counts_by_sample[sample]["expected_retrieved"]),
                          str(counts_by_sample[sample]["retrieved"]),
                          str(counts_by_sample[sample]["detected"]) ])
