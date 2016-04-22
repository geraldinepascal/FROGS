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
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'beta'

import os
import time
import argparse
import subprocess


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def exec_cmd( cmd, output=None ):
    """
    @summary: Execute command line and check status.
    @param cmd: [str] Command line.
    @param output: [str] Path to the file where the stdout will be written.
    """
    if output is None:
        print "\t[Mothur CMD]:\t" + cmd
        subprocess.check_call( cmd, shell=True )
    else:
        print "\t[Mothur CMD]:\t" + cmd + " > " + output
        subprocess.check_call( cmd + " > " + output, shell=True )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == "__main__":
    print "[SOFTWARE]:\tMothur"

    # Manage parameters
    parser = argparse.ArgumentParser(description="Launch mothur workflow.")
    parser.add_argument( '--min-length', type=int, required=True, help='The minimum amplicon length.')
    parser.add_argument( '--max-length', type=int, required=True, help='The maximum amplicon length.')
    parser.add_argument( '--already-contiged', action='store_true', default=False, help='Reads 1 and Reads 2 are already contiged by pair.')
    parser.add_argument( '--preclusters-difference', type=int, default=2, help='The maximum edition distance in pre-clustering.')
    parser.add_argument( '--otu-distance', type=int, default=0.03, help='The OTU maximum distance in clustering.')
    parser.add_argument( '--affi-cutoff', type=int, default=80, help='The bootstrap cutoff in lineage filter.')
    parser.add_argument( '--pcr-start', type=int, help='Start position for amplicon region. This value speedup pipeline by databank restriction.')
    parser.add_argument( '--pcr-end', type=int, help='End position for amplicon region. This value speedup pipeline by databank restriction.')
    parser.add_argument( '--kept-start', type=int, help='In PCR region the start position kept. All sequences must have same size.')
    parser.add_argument( '--kept-end', type=int, help='In PCR region the end position kept. All sequences must have same size.')
    parser.add_argument( '-p', '--nb-cpus', type=int, default=1, help="The maximum number of CPUs used.(default 1)")
    parser.add_argument( '-v', '--version', action='version', version=__version__)
    # Inputs
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('--affiliation-databank-fasta', help='The path to the affiliation databank (format: fasta).')
    group_input.add_argument('--affiliation-databank-tax', help='The path to the affiliation databank (format: mothur tax).')
    group_input.add_argument('--restriction-databank-fasta', required=True, help='The path to sequences used in filters (format: fasta).')
    group_input.add_argument('--restriction-databank-tax', required=True, help='The path to taxonomies used in filters (format: tax).')
    group_input.add_argument('-f', '--input-folder', required=True, help='The path to the input folder.')
    # Outputs
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('--output-biom', required=True, help='The final abundance table with affiliation (format: BIOM).')
    group_output.add_argument('--output-fasta', required=True, help='The final sequences (format: fasta).')
    args = parser.parse_args()


    # Files
    tmp_files = list()
    stability_file = "stability.file"
    args.restriction_databank_fasta = os.path.abspath(args.restriction_databank_fasta)
    args.restriction_databank_tax = os.path.abspath(args.restriction_databank_tax)
    if args.affiliation_databank_fasta is not None:
        args.affiliation_databank_fasta = os.path.abspath(args.affiliation_databank_fasta)
    if args.affiliation_databank_tax is not None:
        args.affiliation_databank_tax = os.path.abspath(args.affiliation_databank_tax)
    args.output_biom = os.path.abspath(args.output_biom)
    args.output_fasta = os.path.abspath(args.output_fasta)

    samples = dict()
    for sample_filename in os.listdir(args.input_folder):
        if sample_filename.endswith(".fastq"):
            sample_name = sample_filename.split(".fastq")[0]
            if sample_name not in samples:
                samples[sample_name] = dict()
            if not args.already_contiged and ("_R2_" in sample_name or "-R2-" in sample_name or sample_name.endswith("_R2") or sample_name.endswith("-R2")):
                samples[sample_name]["R2"] = os.path.join(os.path.abspath(args.input_folder), sample_filename)
            else:
                samples[sample_name]["R1"] = os.path.join(os.path.abspath(args.input_folder), sample_filename)


    #### Move in output directory
    os.chdir(os.path.dirname(args.output_biom))


    #### Write stability file
    reads = list()
    FH_stability = open(stability_file, "w")
    for sample_name in samples:
        sample_line = sample_name + "\t" + samples[sample_name]["R1"]
        reads.append(samples[sample_name]["R1"])
        if not args.already_contiged:
            sample_line += "\t" + samples[sample_name]["R2"]
        else:
            R2 = samples[sample_name]["R1"] + ".RC"
            sample_line += "\t" + R2
            # Create fake R2
            exec_cmd( 'rvc.py ' \
                + '--input ' + samples[sample_name]["R1"] + ' ' \
                + '--output ' + R2 )
            tmp_files.append(R2)
        FH_stability.write(sample_line + "\n")
    FH_stability.close()


    #### Join pairs
    exec_cmd( 'mothur "#make.contigs(' \
        + 'file=' + stability_file + ', ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.fasta
    #   stability.scrap.contigs.fasta
    #   stability.contigs.report
    #   stability.contigs.groups


    #### Start execution time
    start_time = time.time()


    #### Filters
    exec_cmd( 'mothur "#screen.seqs(' \
        + 'fasta=stability.trim.contigs.fasta, ' \
        + 'group=stability.contigs.groups, ' \
        + 'maxambig=0, '\
        + 'maxn=0, ' \
        + 'minlength=' + str(args.min_length) + ', ' \
        + 'maxlength=' + str(args.max_length) + ', ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.fasta
    #   stability.trim.contigs.bad.accnos
    #   stability.contigs.good.groups

    exec_cmd( 'mothur "#unique.seqs(fasta=stability.trim.contigs.good.fasta)"' )
    # Outputs:
    #   stability.trim.contigs.good.names
    #   stability.trim.contigs.good.unique.fasta

    exec_cmd( 'mothur "#count.seqs(' \
        + 'name=stability.trim.contigs.good.names, ' \
        + 'group=stability.contigs.good.groups)"' )
    # Outputs:
    #   stability.trim.contigs.good.count_table


    #### Align on targeted region
    # Restrict silva to primer area (V3V4: end 25319; V4V4: end 25319)
    exec_cmd( 'cp ' + args.restriction_databank_fasta + ' restriction_db.fasta' )
    exec_cmd( 'mothur "#pcr.seqs(' \
        + 'fasta=restriction_db.fasta, ' \
        + 'keepprimer=T, ' \
        + 'start=' + str(args.pcr_start) + ', ' \
        + 'end=' + str(args.pcr_end) + ', ' \
        + 'keepdots=F, ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   silva.bacteria.pcr.fasta

    # Align sequences on databank
    exec_cmd( 'mothur "#align.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.fasta, ' \
        + 'reference=restriction_db.pcr.fasta, ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.align
    #   stability.trim.contigs.good.unique.align.report
    exec_cmd( 'mothur "#summary.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.align, ' \
        + 'count=stability.trim.contigs.good.count_table, ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.summary

    # Remove sequence unaligned in right area
    exec_cmd( 'mothur "#screen.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.align, ' \
        + 'count=stability.trim.contigs.good.count_table, ' \
        + 'summary=stability.trim.contigs.good.unique.summary, ' \
        + 'start=' + str(args.kept_start) + ', ' \
        + 'end=' + str(args.kept_end) + ', ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.summary
    #   stability.trim.contigs.good.unique.good.align
    #   stability.trim.contigs.good.unique.bad.accnos
    #   stability.trim.contigs.good.good.count_table
    exec_cmd( 'mothur "#filter.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.good.align, ' \
        + 'vertical=T, ' \
        + 'trump=., ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.filter
    #   stability.trim.contigs.good.unique.good.filter.fasta

    exec_cmd( 'mothur "#unique.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.fasta, ' \
        + 'count=stability.trim.contigs.good.good.count_table)"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.count_table
    #   stability.trim.contigs.good.unique.good.filter.unique.fasta


    #### Pre-cluster
    exec_cmd( 'mothur "#pre.cluster(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, ' \
        + 'count=stability.trim.contigs.good.unique.good.filter.count_table, ' \
        + 'diffs=' + str(args.preclusters_difference) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.mock_v4.map


    #### Remove singleton


    #### Remove chimera
    exec_cmd( 'mothur "#chimera.uchime(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, ' \
        + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, ' \
        + 'dereplicate=t, ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
    if os.path.exists("stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos"): # mothur version < 1.35.0
        exec_cmd( 'ln -sf ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table' )
        exec_cmd( 'ln -sf ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.chimeras ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras' )
        exec_cmd( 'ln -sf ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos ' \
            + 'stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos' )
    exec_cmd( 'mothur "#remove.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, ' \
        + 'accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta


    #### Affiliation seq and reduce lineage
    exec_cmd( 'cp ' + args.restriction_databank_tax + ' restriction_db.tax' )
    classif_lineage_time = time.time()
    exec_cmd( 'mothur "#classify.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, ' \
        + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, ' \
        + 'reference=restriction_db.fasta, ' \
        + 'taxonomy=restriction_db.tax, ' \
        + 'cutoff=' + str(args.affi_cutoff) + ')"' )################################### + 'processors=' + str(args.nb_cpus) + ')"'
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.restriction_db.wang.taxonomy
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.restriction_db.wang.tax.summary
    classif_lineage_time = time.time() - classif_lineage_time
    exec_cmd( 'mothur "#remove.lineage(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, ' \
        + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, ' \
        + 'taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.restriction_db.wang.taxonomy, ' \
        + 'taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.restriction_db.wang.pick.taxonomy
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table


    #### Clustering
    exec_cmd( 'mothur "#dist.seqs(' \
        + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, ' \
        + 'cutoff=0.20, ' \
        + 'processors=' + str(args.nb_cpus) + ')"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist
    exec_cmd( 'mothur "#cluster(' \
        + 'column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, ' \
        + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table)"' )
    # Outputs:
    #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list


    if args.affiliation_databank_fasta is None or args.affiliation_databank_tax is None:
        #### BIOM
        exec_cmd( 'mothur "#make.shared(' \
            + 'list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, ' \
            + 'label=' + str(args.otu_distance) + ')"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.mock_v4.rabund

        exec_cmd( 'mothur "#make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared)"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom

        exec_cmd( 'ln -sf stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom ' + args.output_biom )

        exec_cmd( 'mothur "#get.oturep(' \
            + 'method=abundance, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, ' \
            + 'list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, ' \
            + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, ' \
            + 'cutoff=' + str(args.otu_distance) + ')"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.count_table
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta

        # Degap and rename OTU
        exec_cmd( 'mothurDeGapSeeds.py ' \
            + '--input stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.' + str(args.otu_distance) + '.rep.fasta ' \
            + '--output stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.' + str(args.otu_distance) + '.rep.degap.fasta' )

        # Add reference ID in seeds descriptions
        exec_cmd( 'mothurAddSeedRef.py ' \
            + "--input stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list." + str(args.otu_distance) + ".rep.degap.fasta "  \
            + "--reads " + " ".join(reads) + " " \
            + "--trimmed-reads stability.trim.contigs.good.unique.good.filter.fasta " \
            + "--output " + args.output_fasta )
    else:
        exec_cmd( 'cp ' + args.affiliation_databank_fasta + ' affiliation_db.fasta' )
        exec_cmd( 'cp ' + args.affiliation_databank_tax + ' affiliation_db.tax' )

        #### Affiliation
        exec_cmd( 'mothur "#classify.seqs(' \
            + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, ' \
            + 'reference=affiliation_db.fasta, ' \
            + 'taxonomy=affiliation_db.tax, ' \
            + 'cutoff=0)"' )################################### + 'processors=' + str(args.nb_cpus) + ')"'
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.affiliation_db.wang.taxonomy
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.affiliation_db.wang.tax.summary
        exec_cmd( 'mothur "#classify.otu('
            + 'list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, ' \
            + 'taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.affiliation_db.wang.taxonomy, ' \
            + 'label=' + str(args.otu_distance) + ')"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.tax.summary


        #### BIOM
        exec_cmd( 'mothur "#make.shared(' \
            + 'list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, ' \
            + 'label=' + str(args.otu_distance) + ')"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.*.rabund

        exec_cmd( 'mothur "#make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, ' \
            + 'constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.' + str(args.otu_distance) + '.cons.taxonomy)"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom

        exec_cmd( 'ln -sf stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.biom ' + args.output_biom )

        exec_cmd( 'mothur "#get.oturep(' \
            + 'method=abundance, ' \
            + 'count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, ' \
            + 'list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, ' \
            + 'fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, ' \
            + 'cutoff=' + str(args.otu_distance) + ')"' )
        # Outputs:
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.count_table
        #   stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.rep.fasta

        # Degap and rename OTU
        exec_cmd( 'mothurDeGapSeeds.py ' \
            + '--input stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.' + str(args.otu_distance) + '.rep.fasta ' \
            + '--output stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.' + str(args.otu_distance) + '.rep.degap.fasta' )

        # Add reference ID in seeds descriptions
        exec_cmd( 'mothurAddSeedRef.py ' \
            + "--input stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list." + str(args.otu_distance) + ".rep.degap.fasta "  \
            + "--reads " + " ".join(reads) + " " \
            + "--trimmed-reads stability.trim.contigs.good.unique.good.filter.fasta " \
            + "--output " + args.output_fasta )


    #### Print execution time
    end_time = time.time()
    print "\t[Mothur EXEC_TIME]:\t" + str((end_time - start_time) - classif_lineage_time)


    #### Cleanning tmp
    exec_cmd( 'rm -f mothur.*.logfile ' + " ".join(tmp_files) )
