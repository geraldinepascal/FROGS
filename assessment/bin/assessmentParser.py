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

__author__ = 'Plateforme bioinformatique Toulouse - Sigenae Jouy en Josas'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import os
import json
import argparse
import subprocess


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def parseOTUMetrics(OTU_results_path):
    OTU_metrics = dict()
    if os.path.exists(OTU_results_path):
        FH_results = open(OTU_results_path)
        is_samples_section = False
        titles = None
        for line in FH_results:
            line = line.strip()
            if line.startswith("#Sample"): # '#Sample    After_simu    Distincts_after_simu    Expected_retrieved    Retrieved    Detected'
                is_samples_section = True
                titles = line.split()
            elif is_samples_section:
                line_fields = line.split()
                OTU_metrics[line_fields[0]] = dict()
                for title_idx, val in enumerate(line_fields[1:]):
                    OTU_metrics[line_fields[0]][titles[title_idx+1]] = int(val)
        FH_results.close()
    return OTU_metrics


def parseAffiliationsMetrics(affi_results_path):
    """
sample01-20sp-Powerlaw
#Rank    Divergence (%)    Common    Real specific    Frogs specific
Domain    0.0    1    0    0
Phylum    0.351903607876    5    0    0
Class    0.424546670058    8    0    0
Order    0.431419145727    11    0    0
Family    0.560344622176    15    0    0
Genus    0.634757175613    18    0    0
Species    0.73745107518    20    0    3

sample02-20sp-Powerlaw
    """
    affi_metrics = dict()
    if os.path.exists(affi_results_path):
        FH_results = open(affi_results_path)
        is_new_section = True
        sample = None
        for line in FH_results:
            line = line.strip()
            if line == "":
                is_new_section = True
            elif is_new_section:
                is_new_section = False
                sample = line
                affi_metrics[sample] = dict()
            elif not line.startswith("#") and not line.startswith("DIFF") and not line.startswith("GRINDER"):
                line_fields = line.split("\t")
                affi_metrics[sample][line_fields[0]] = float(line_fields[1])
        FH_results.close()
    return affi_metrics


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
    parser.add_argument( '--pipelines', type=str, nargs='+', default=["frogs", "uparse", "mothur", "qiime"], help='launched pipelines. Example: with "--pipelines frogs" only the frogs is launched. [Default: frogs, uparse, mothur].' )
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
                if not pipeline == "frogs":
                    workflows_types.append(pipeline+"-FROGS-affi")
            for current_workflow in workflows_types:
                workflows_metics = {"clusters": dict(), "affiliations": dict()}
                for current_nb_sp in args.nb_sp:
                    workflows_metics["clusters"][current_nb_sp] = dict()
                    workflows_metics["affiliations"][current_nb_sp] = dict()
                    for dataset_idx in args.datasets:
                        dataset_in_dir = os.sep.join([args.assessment_directory, str(current_nb_sp) + "sp", "dataset_" + str(dataset_idx), current_primers, current_distribution])
                        wf_directory = os.path.join(dataset_in_dir, current_workflow.split("-")[0])
                        clstr_results_path = os.path.join(wf_directory, current_workflow + "_OTUResults.txt")
                        affi_results_path = os.path.join(wf_directory, current_workflow + "_affiResults.txt")
                        workflows_metics["clusters"][current_nb_sp]["datatset_" + str(dataset_idx)] = parseOTUMetrics(clstr_results_path)
                        workflows_metics["affiliations"][current_nb_sp]["datatset_" + str(dataset_idx)] = parseAffiliationsMetrics(affi_results_path)
                print current_workflow + "_" + current_primers + "_" + current_distribution + " processed"
                wf_results_path = os.path.join(args.output_directory, current_workflow + "_" + current_primers + "_" + current_distribution + ".json")
                FH_wf_results = open(wf_results_path, "w")
                FH_wf_results.write( json.dumps(workflows_metics, sort_keys=True, indent=3) )
                FH_wf_results.close()
