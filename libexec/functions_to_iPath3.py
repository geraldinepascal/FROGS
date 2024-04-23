import sys
import os
from collections import defaultdict
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests
from IPython.display import SVG
import svgutils.compose as sc

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib/external-lib"))
PATHWAY_FILE= os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "frogsfunc_suppdata/Kegg_pathways_to_functions.tsv"))
sys.path.append(LIB_DIR)

import ipath

def generate_ipath_input( file ):
    out_str = str()
    with open(file) as FH_in:
        for li in FH_in:
            if li.startswith("ASV"):
                continue
            else:
                li = li.strip().split(' ')
                out_str += ' '.join([ li[0], li[1], li[2] ]) + '\n'

    return out_str


def dict_init_path_func( file ):
    path_to_count = dict()
    func_to_paths = defaultdict(list)
    tot_func = 0
    FH_in = open(file).readlines()
    for li in FH_in:
        li = li.strip().split('\t')
        tot_func += int(li[2])
        for func in li[3].split(';'):
            func_to_paths[func].append(li[0])
    for li in FH_in:
        li = li.strip().split('\t')
        # init map id with map name, nb func, tot nb func, 0 and 0 for init diff func retrieved
        path_to_count[li[0]] = [li[1], int(li[2]), tot_func, 0, 0]

    return path_to_count, func_to_paths


def retrieve_diff_functions( file, dict_path, dict_func ):
    tot_func = 0
    with open(file) as FH_in:
        for li in FH_in:
            li = li.strip().split('\t')
            if li[0] in dict_func:
                for path in dict_func[li[0]]:
                    tot_func += 1
                    if path in dict_path:
                        dict_path[path][3] += 1
    for path in dict_path:
        dict_path[path][4] = tot_func
    
    return dict_path

def process_stats( results, p_thresh ):

    to_adj = list()
    sign_results = dict()
    cur = 0

    for path_id, infos in results.items():
        ratio, p = stats.fisher_exact([[infos[2], infos[1]], [infos[4], infos[3]]])

        if p < p_thresh:
            if ratio > 1:
                sens = "Over_represented"
            else:
                sens = "Under_represented"
            sign_results[cur] = [infos[0], infos[2], infos[1], infos[4], infos[3], sens, ratio, p ]
            to_adj.append(p)
            cur += 1
    
    p_adjusted = multipletests(to_adj, method='bonferroni', alpha=p_thresh)

    for j in range(len(list(p_adjusted[1]))):
        padj = list(p_adjusted[1])[j]
        if padj < p_thresh:
            sign_results[j].append(padj)

            if padj < 0.001 :
                sign = "***"
            elif 0.001 < padj < 0.01:
                sign = "**"
            elif 0.01 < padj < 0.05:
                sign = "*"
            elif 0.05 < padj < 0.1:
                sign = "."
            else:
                sign = " "
            sign_results[j].append(sign)

        else:
            del sign_results[j]

    return sign_results

def write_output(sign_results, output):

	FH_out = open(output, 'wt')

	FH_out.write('\t'.join(["Pathway_name", "Nb_total_functions_all_pathways",\
             "Nb_total_functions_in_pathway", "Nb_total_diff_functions",\
             "Nb_diff_functions_in_pathway", "Sens_of_the_enrichment", "Ratio", "p-value", "adjusted_p-value", "sign"]) + "\n" )

	for id, signs in sign_results.items():
		FH_out.write("\t".join(map(str,signs))+"\n")

up_file = "over.tsv"
down_file = "under.tsv"
# p_thresh = 0.05

# path_to_count, func_to_paths = dict_init_path_func(PATHWAY_FILE)

# results_path_to_count = retrieve_diff_functions( up_file, path_to_count, func_to_paths )

# sign_results = process_stats(results_path_to_count, p_thresh)

# write_output( sign_results, "sign_path.tsv")

for file in [up_file, down_file]:
    str_out = str()
    str_out = generate_ipath_input( file)
    ipath.get_map(str_out, os.path.basename(file).split('.')[0])
