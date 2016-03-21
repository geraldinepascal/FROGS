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
__status__ = 'prod'

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
    parser = argparse.ArgumentParser( description='De gap sequence and replace sequence ID by OTU ID.' )
    parser.add_argument( '-v', '--version', action='version', version=__version__ )
    parser.add_argument( '-i', '--input', required=True, help='Sequences file from mothur get.oturep (format: fasta).' ) 
    parser.add_argument( '-o', '--output', required=True, help='Output file (format: fasta).' )
    args = parser.parse_args()
    
    # Process
    FH_in = FastaIO(args.input)
    FH_out = FastaIO(args.output, "w")
    for record in FH_in:
        seed_id = record.id
        record.id = record.description.split("|")[0]
        record.description += " seed_id=" + seed_id
        record.string = record.string.replace("-", "").replace(".", "")
        FH_out.write(record)
    FH_in.close()
    FH_out.close()
