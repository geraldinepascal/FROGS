#!/usr/bin/env python2.7
#
# Copyright (C) 2014 INRA
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

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1.1'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'prod'


import os
import sys
import json
import shutil
import urllib
import tarfile
import argparse
# Galaxy dependencies
import galaxy.model # need to import model before sniff to resolve a circular import dependency
from galaxy.datatypes import sniff
from galaxy.datatypes.registry import Registry
from galaxy import util


##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_metadata(path, dataset):
    """
    @summary: Writes galaxy metadata for the dataset.
    @param path: [str] Path to the metadata file.
    @param dataset: [Bunch] The dataset to process.
    """
    FH = open( path, 'w' )
    info = dict( type = 'dataset',
                 dataset_id = dataset.dataset_id,
                 ext = "data",
                 stdout = "uploaded archive",
                 name = dataset.name )
    if dataset.get('uuid', None) is not None:
        info['uuid'] = dataset.get('uuid')
    FH.write( json.dumps( info ) + "\n" )
    FH.close()

def url_upload(dataset):
    """
    @summary: Upload an URL file and update dataset information.
    @param dataset: [Bunch] The dataset for the uploaded file.
    """
    url = dataset.path
    page = urllib.urlopen( url )
    temp_name, dataset.is_multi_byte = sniff.stream_to_file( page, prefix='url_paste', source_encoding=util.get_charset_from_http_headers( page.headers ) )
    dataset.path = temp_name
    dataset.name = url.split('/')[-1]

def safe_dict(d):
    """
    Recursively clone json structure with UTF-8 dictionary keys
    http://mellowmachines.com/blog/2009/06/exploding-dictionary-with-unicode-keys-as-python-arguments/
    """
    if isinstance(d, dict):
        return dict([(k.encode('utf-8'), safe_dict(v)) for k,v in d.iteritems()])
    elif isinstance(d, list):
        return [safe_dict(x) for x in d]
    else:
        return d

class OutputParameter(argparse.Action):
    """
    @summary : Argparse parameter for output.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        # Set parser
        output = getattr(namespace, self.dest)
        if output is None:
            output = dict()
        # Retrieve params
        for option in values:
            id, filepath = option.split(":")
            output[id] = filepath
        setattr(namespace, self.dest, output)

def process(param_file, output_paths):
    """
    @summary: Upload an archive in galaxy dataset.
    @param param_file: [str] Path to the galaxy paramfile.
    @param output_paths: [dict] The output paths by dataset ID.
    """
    try:
        for line in open( param_file, 'r' ):
            dataset = json.loads( line )
            dataset = util.bunch.Bunch( **safe_dict( dataset ) )
            if dataset.type == 'url':
                url_upload( dataset )
            # Check file
            if not tarfile.is_tarfile(dataset.path):
                os.remove( dataset.path )
                raise Exception("The archive '" + dataset.name + "' is not a tar file.")
            # Move file
            output_path = output_paths[str(dataset.dataset_id)]
            shutil.move( dataset.path, output_path )
            dataset.path = output_path
            # Write metadata
            write_metadata( 'galaxy.json', dataset )
    finally:
        os.remove( param_file )


##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################
if __name__ == '__main__':
    # Manage parameters
    parser = argparse.ArgumentParser(
        description='Upload an archive file (.tar or .tar.gz).'
    )
    parser.add_argument( '-p', '--param-file', type=str, required=True, help='The parameter file from UploadToolAction.' )
    parser.add_argument( '-o', '--output', type=str, action=OutputParameter, nargs="+", metavar=("DATASET_ID:OUTPUT_DATASET_PATH", ""), help='The output path for each element in paramfile.' )
    args = parser.parse_args()

    # Process
    process( args.param_file, args.output )