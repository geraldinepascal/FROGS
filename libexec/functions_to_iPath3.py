import sys
import os
from IPython.display import SVG
import svgutils.compose as sc

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "../lib/external-lib"))
sys.path.append(LIB_DIR)

import ipath

def generate_ipath_input( file ):
    out_str = str()
    with open(file) as FH_in:
        for li in FH_in:
            if li.startswith("OTU"):
                continue
            else:
                li = li.strip().split('\t')
                out_str += ' '.join([ li[0], li[1], li[2] ]) + '\n'
    return out_str

up_file = "../../tools/deseq2_visualisation/over.tsv"
down_file = "../../tools/deseq2_visualisation/under.tsv"

for file in [up_file, down_file]:
    str_out = str()
    str_out = generate_ipath_input( file)
    ipath.get_map(str_out, os.path.basename(file).split('.')[0])
