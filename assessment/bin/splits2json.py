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

__author__ = 'Plateforme bioinformatique Toulouse - Sigenae'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'beta'

import os
import json
import argparse
import subprocess


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def parseOTUSplits(OTU_results_path):
    nb_splits = None
    if os.path.exists(OTU_results_path):
        FH_results = open(OTU_results_path)
        is_dataset_section = False
        titles = None
        for line in FH_results:
            line = line.strip()
            if line.startswith("#After_simu"): # '#After_simu    Dictincts_after_simu    Expected_retrieved    Retrieved    Detected    Splits'
                is_dataset_section = True
                titles = line.split()
            elif is_dataset_section and line == "":
                is_dataset_section = False
            elif is_dataset_section:
                line_fields = line.split()
                splits_idx = titles.index("Splits")
                nb_splits = int(line_fields[splits_idx])
        FH_results.close()
    return nb_splits


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    # Manage parameters
    parser = argparse.ArgumentParser(description="Format assessment results in JSON files.")
    parser.add_argument( '--nb-sp', type=int, nargs='+', default=[20, 100, 200, 500, 1000], help='Numbers of species used in assessment. Example: with "--nb-sp 20 100" only the results for 20sp and 100sp are converted. [Default: 20, 100, 200, 500, 1000]' )
    parser.add_argument( '--datasets', type=int, nargs='+', default=[1, 2, 3, 4, 5], help='IDs of datasets used in assessment. Example: with "--datasets 1 5" only the results for dataset_1 and dataset_5 are converted. [Default: 1, 2, 3, 4, 5]' )
    parser.add_argument( '--primers', type=str, nargs='+', default=["V4V4", "V3V4"], help='Primers used in assessment. Example: with "--primers V4V4" only the results for V4V4 are converted. [Default: V4V4, V3V4]' )
    parser.add_argument( '--distribution-laws', type=str, nargs='+', default=["powerlaw", "uniform"], help='Distribution laws used in assessment. Example: with "--distribution-laws powerlaw" only the results for powerlaw are converted. [Default: powerlaw, uniform].' )
    parser.add_argument( '--pipelines', type=str, nargs='+', default=["frogs", "uparse", "mothur", "qiime"], help='launched pipelines. Example: with "--pipelines frogs" only the frogs is launched. [Default: frogs, uparse, mothur]. For qiime do not forget to configure your ~/.qiime_config file with your reference database' )
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument( '--assessment-directory', required=True, help='Path to the output directory of the assessment.' )
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument( '--output-directory', required=True, help='Path to the output directory.' )
    args = parser.parse_args()

    # Process
    for current_primers in args.primers:
        for current_distribution in args.distribution_laws:
            workflows_types = list()
            for pipeline in args.pipelines : 
                workflows_types.append(pipeline)
                workflows_types.append(pipeline+"-FROGS-affi")
            for current_workflow in workflows_types:
                workflows_metics = dict()
                for current_nb_sp in args.nb_sp:
                    workflows_metics[current_nb_sp] = dict()
                    for dataset_idx in args.datasets:
                        dataset_in_dir = os.sep.join([args.assessment_directory, str(current_nb_sp) + "sp", "dataset_" + str(dataset_idx), current_primers, current_distribution])
                        wf_directory = os.path.join(dataset_in_dir, current_workflow.split("-")[0])
                        clstr_results_path = os.path.join(wf_directory, current_workflow + "_OTUResults.txt")
                        workflows_metics[current_nb_sp]["datatset_" + str(dataset_idx)] = parseOTUSplits(clstr_results_path)
                print current_workflow + "_" + current_primers + "_" + current_distribution + " processed"
                wf_results_path = os.path.join(args.output_directory, current_workflow + "_" + current_primers + "_" + current_distribution + "_splits.json")
                FH_wf_results = open(wf_results_path, "w")
                FH_wf_results.write( json.dumps(workflows_metics, sort_keys=True, indent=3) )
                FH_wf_results.close()
