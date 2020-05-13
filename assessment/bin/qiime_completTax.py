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

__author__ = 'Sigenae INRA Jouy en Josas'
__copyright__ = 'Copyright (C) 2016 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'prod'

import argparse

def completTax(in_tax,out_tax):
    utax_ranks = ["d:", "p:", "c:", "o:", "f:", "g:"]
    max_rank=0
    tmp_lines=list()
    FH_in = open(in_tax)
    FH_out = open(out_tax,"w")
    for line in FH_in:
        taxonomy=line.split("\t")[1] 
        new_tax = taxonomy if not taxonomy.endswith(";") else taxonomy[:-1]
        # incomplet taxonomy
        if not taxonomy == "Unassigned"  and taxonomy.startswith("d:") : #utax  database
            max_rank=6
            if len(new_tax.split(";")) < 6:
                last_rank = taxonomy.split(";")[-1][0:2]
                idx_rank = utax_ranks.index(last_rank)
                for i in xrange(idx_rank+1, len(utax_ranks)):
                    new_tax += ";"+utax_ranks[i]+"unknown_taxa"
        elif not taxonomy == "Unassigned" : # silva database
            max_rank=7
            idx_rank =  len(taxonomy.split(";"))
            for i in xrange(idx_rank, 7):
                new_tax += ";"+"unknown_taxa"
        elif taxonomy=="Unassigned":
            if max_rank >0:
                if tmp_lines != []:
                    for l in tmp_lines : 
                        tmp_taxonomy = line.split("\t")[1] 
                        tmp_new_tax = tmp_taxonomy+";"
                        tmp_new_tax *= max_rank
                        tmp_new_tax = tmp_new_tax[:-1]
                        l = l.replace(tmp_taxonomy,tmp_new_tax)
                        FH_out.write(l)   
                    tmp_lines = []                    
                new_tax += ";"
                new_tax *= max_rank
                new_tax = new_tax[:-1]
            else:
                tmp_lines.append(line)
                continue
        line=line.replace(taxonomy,new_tax)
        FH_out.write(line)
    if tmp_lines != []:
        for l in tmp_lines : 
            tmp_taxonomy = line.split("\t")[1] 
            tmp_new_tax = tmp_taxonomy+";"
            tmp_new_tax *= max_rank
            tmp_new_tax = tmp_new_tax[:-1]
            l = l.replace(tmp_taxonomy,tmp_new_tax)
            FH_out.write(l)   
    FH_in.close()
    FH_out.close()
    
if __name__ == "__main__":

    # Manage parameters
    parser = argparse.ArgumentParser(description="Complete qiime taxonomy to have always 6 ranks (for utax) or 7 ranks (for silva).")
    group_input = parser.add_argument_group('Inputs')
    group_input.add_argument('-i', '--input-tax', required=True, help='The Qiime output taxonomy file')
    group_output = parser.add_argument_group('Outputs')
    group_output.add_argument('-o', '--output-tax', required=True, help='The completed taxonomy output file')
    args = parser.parse_args()
    
    completTax(args.input_tax,args.output_tax)
