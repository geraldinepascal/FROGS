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

from unicodedata import name
import mechanize
import argparse
import requests
import gzip
import re
from bs4 import BeautifulSoup

### URL PATTERN DATABASES
JGI_URL = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid='
ITS_URL = 'https://mycocosm.jgi.doe.gov/###ID/###ID.home.html'
TAXO_BROWSER = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
###

### Classes 

###

### Functions
# tools
def make_a_soup(url):
    html_text = requests.get(url).text
    soup = BeautifulSoup(html_text, 'html.parser')
    return soup

def run_form(to_search, name):
    to_search.open(TAXO_BROWSER)
    to_search.select_form(nr=0)
    to_search.form['srchmode'] = ['1']
    to_search['name'] = name
    req = to_search.submit()
    for li in req:
        li = str(li)
        if "Did you mean" in li:
            name = li.split('alt="score=0">')[1].split('</a>')[0]
            return run_form(to_search, name)
    return name

def is_gzipped(path):
    return path.endswith(".gz")

# parsing input files
def read_alignment_file(alignment_fi):
    opener = gzip.open if is_gzipped(alignment_fi) else open
    with opener(alignment_fi, "rt") as fi:
        for li in fi:
            li = li.strip()
            if li.startswith('>'):
                yield li.split('>')[-1].split('-cluster')[0]

# Parsing html pages
def parse_jgi_html(id, jgi_url, output_file):
    ###
    url = jgi_url + id
    soup = make_a_soup(url)
    # Species name parser
    try:
        name = soup.find(class_="subhead", text="Organism Name")\
            .find_next_sibling('td').text
    except:
        try:
            removed = soup.find(id="content_other")
            id = str(removed.b).split('taxon_oid=')[1].split('">')[0]
            return parse_jgi_html(id, jgi_url, output_file)
        except:
            name = "unknown"

    # lineage parser
    try:
        # gtdb-tk clasification
        lineage = ";".join([cur_rank.split('__')[-1] for cur_rank in \
            soup.find(class_="subhead", text="GTDB-tk Lineage").find_next_sibling('td').text.split(';')])
    except:
        # else, retrieve default classification
        try:
            lineage = ";".join([cur_rank.get_text() for cur_rank in  \
                soup.find(class_="subhead", text="Lineage").find_next_sibling('td').find_all('a')])
        except:
            lineage = ",,,,,,"
            
    output_file.write(f'{id}\t{name}\t{lineage}\n')

def parse_its_html(id, jgi_url, to_search, output_file):
    ###
    url = jgi_url.replace('###ID',id)
    soup = make_a_soup(url)
    REGEX_VERSION = "v[0-9]{1,2}(\.[0-9]{1,2})?"
    # Species name parser
    name = soup.find(class_="organismName").text.replace("Home • ","")
    if re.search(REGEX_VERSION, name):
        name = " ".join(name.split()[:-1])
                
def parse_arguments():
    # Manage parameters.
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--ref_alignment", required = True, \
     help = "PICRUSt2 reference alignment file")

    parser.add_argument("-c", "--category", required = True, \
     help = "Category of database to be searched. 16S : JGI databbase, ITS : JGI Mycocosm database" , choices=["16S", "ITS"])

    parser.add_argument('-o', '--output_file', \
    default="JGI_ID_to_taxonomy.txt", help = 'Output file with JGI metadata.')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    ids = read_alignment_file(args.ref_alignment)

    if args.category == "ITS":
            to_search = mechanize.Browser()

    output_fi = open(args.output_file,'w')
    for id in ids:
        if args.category == "16S":
            parse_jgi_html(id, JGI_URL, output_fi)
        elif args.category == "ITS":
            parse_its_html(id, ITS_URL, to_search, output_fi)
    output_fi.close()

if __name__ == '__main__':
    main()