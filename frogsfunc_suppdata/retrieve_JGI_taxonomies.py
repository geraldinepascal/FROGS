#!/usr/bin/env python3
#
# Copyright (C) 2022 INRA
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
__author__ = 'Vincent Darbot INRAE - GENPHYSE'
__copyright__ = 'Copyright (C) 2022 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'dev'

import argparse
import sys
import re
import os
import requests
from bs4 import BeautifulSoup

JGI_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid='

def read_alignment_file(alignment_fi):
    with open(alignment_fi) as fi:
        for li in fi:
            li = li.strip()
            if li.startswith('>'):
                yield li.split('>')[-1].split('-cluster')[0]

def parse_jgi_html(id, jgi_url, output_file):
    ###
    url = jgi_url + id
    html_text = requests.get(url).text
    soup = BeautifulSoup(html_text, 'html.parser')
    # Species name parser
    try:
        name = soup.find(class_="subhead", text="Organism Name")\
            .find_next_sibling('td').text
    except:
        name = "unknown"
    # lineage parser
    try:
        # gtdb-tk clasification
        lineage = ";".join([cur_rank.split('__')[-1] for cur_rank in \
            soup.find(class_="subhead", text="GTDB-tk Lineage").find_next_sibling('td').text.split(';')])
    except NoLineageException:
        # else, retrieve default classification
        lineage = ";".join([cur_rank.get_text() for cur_rank in  \
            soup.find(class_="subhead", text="Lineage").find_next_sibling('td').find_all('a')])
    except:
        lineage = ";;;;;;"
    ###
    output_file.write(f'{id}\t{name}\t{lineage}\n')

def parse_arguments():
    # Manage parameters.
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--ref_alignment", required = True, \
     help = "PICRUSt2 reference alignment file")

    parser.add_argument('-o', '--output_file', \
    default="JGI_ID_to_taxonomy.txt", help = 'Output file with JGI metadata.')

    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    ids = read_alignment_file(args.ref_alignment)
    
    output_fi = open(args.output_file,'w')
    for id in ids:
        parse_jgi_html(id, JGI_URL, output_fi)
    output_fi.close()


if __name__ == '__main__':
    main()