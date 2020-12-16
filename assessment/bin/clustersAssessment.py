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
__version__ = '1.2.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'beta'

import re
import os
import argparse
from frogsSequenceIO import *
from frogsBiom import BiomIO


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
    uniq_id = dict()
    uniq_id_by_sample = dict()

    # Get sequences
    ids_by_sequences = dict()
    fh_sequences = FastaIO( fasta )
    for record in fh_sequences:
        if record.string not in ids_by_sequences:
            ids_by_sequences[record.string] = list()
        ids_by_sequences[record.string].append( record.id )
    fh_sequences.close()

    # Get uniq in all dataset
    for sequence in ids_by_sequences:
        for reference_id in ids_by_sequences[sequence]:
            uniq_id[reference_id] = ids_by_sequences[sequence][0]

    # Get uniq by sample
    for sample_name in references_by_sample:
        uniq_id_by_sample[sample_name] = dict()
        for reference_id in references_by_sample[sample_name]:
            uniq_id_by_sample[sample_name][reference_id] = uniq_id[reference_id]

    return uniq_id, uniq_id_by_sample


def get_retrieved_in_dataset(reference_by_obs_id, references, uniq_id):
    nb_detected = 0
    retrieved = dict()
    expected_retrieved = dict()
    nb_splits = 0
    for obs_id in reference_by_obs_id:
        nb_detected += 1
        if not "," in reference_by_obs_id[obs_id]: # Is not a chimera
            ref_id = reference_by_obs_id[obs_id]
            if ref_id in retrieved:
                nb_splits += 1
            else:
                retrieved[ref_id] = 1
                if ref_id in references:
                    expected_retrieved[ref_id] = 1
    # Uniq sequence for retrieved
    uniq_retrieved = set()
    for ref_id in retrieved:
        uniq_retrieved.add( uniq_id[ref_id] )
    # Add to split the references with the same sequence and retrieved separately
    nb_splits += len(retrieved.keys()) - len(uniq_retrieved)
    # Uniq sequence for retrieved
    uniq_expected_retrieved = set()
    for ref_id in expected_retrieved:
        uniq_expected_retrieved.add( uniq_id[ref_id] )
    # Results
    return {
        "expected_retrieved": len(uniq_expected_retrieved),
        "retrieved": len(uniq_retrieved),
        "detected": nb_detected,
        "splits": nb_splits
    }


def get_retrieved_by_sample( biom_file, reference_by_obs_id, references_by_sample, uniq_id, uniq_id_by_sample ):
    counts_by_sample = dict()
    biom = BiomIO.from_json( biom_file )
    for sample_name in biom.get_samples_names():
        nb_detected = 0
        retrieved = dict()
        expected_retrieved = dict()
        for obs in biom.get_observations_by_sample( sample_name ):
            nb_detected += 1
            if not "," in reference_by_obs_id[obs['id']]: # Is not a chimera
                ref_id = reference_by_obs_id[obs['id']]
                retrieved[ref_id] = 1
                if ref_id in references_by_sample[sample_name]:
                    expected_retrieved[ref_id] = 1
        # Uniq sequence for retrieved
        uniq_retrieved = set()
        for ref_id in retrieved:
            uniq_retrieved.add( uniq_id[ref_id] )
        # Uniq sequence for retrieved
        uniq_expected_retrieved = set()
        for ref_id in expected_retrieved:
            uniq_expected_retrieved.add( uniq_id_by_sample[sample_name][ref_id] )
        # Results
        counts_by_sample[sample_name] = {
            "detected": nb_detected,
            "retrieved": len(uniq_retrieved),
            "expected_retrieved": len(uniq_expected_retrieved)
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
    uniq_id, uniq_id_by_sample = get_uniq(args.origin, references_by_sample)

    # Get OTU in results
    reference_by_obs_id = dict()
    fh_sequences = FastaIO( args.fasta )
    for record in fh_sequences:
        if ";size=" in record.id:
            record.id = record.id.split(";size=", 1)[0]
        reference_by_obs_id[record.id] = re.search("reference=([^\s]+)", record.description).group(1)
    fh_sequences.close()
    counts = get_retrieved_in_dataset( reference_by_obs_id, references, uniq_id )
    counts_by_sample = get_retrieved_by_sample( args.biom, reference_by_obs_id, references_by_sample, uniq_id, uniq_id_by_sample )

    # Output
    print "#After_simu\tDictincts_after_simu\tExpected_retrieved\tRetrieved\tDetected\tSplits"
    print "\t".join([ str(len(references)),
                      str(len(set(uniq_id.values()))),
                      str(counts["expected_retrieved"]),
                      str(counts["retrieved"]),
                      str(counts["detected"]),
                      str(counts["splits"]) ])
    print ""
    
    print "#Sample\tAfter_simu\tDistincts_after_simu\tExpected_retrieved\tRetrieved\tDetected"
    for sample in samples:
        print "\t".join([ sample, 
                          str(len(references_by_sample[sample])),
                          str(len(set(uniq_id_by_sample[sample].values()))),
                          str(counts_by_sample[sample]["expected_retrieved"]),
                          str(counts_by_sample[sample]["retrieved"]),
                          str(counts_by_sample[sample]["detected"]) ])
