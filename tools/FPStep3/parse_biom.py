#Script pour parser le fichier biom
"""
Format souhaiter:
id  sample1 sample2
cluster_1   sum_obs1 sum_obs2
clsuter_2   sum_obs1    sum_obs2
"""
#Le biom de FROGS sort comme un json
#Le biom que prend Picrust2 est sous forme de tsv
#intégré FROGS_BIOM_To_TSV  récupérer le TSV et le donner à FPSTep3
import os
import sys
import json
import re
def test(biom, out_file1):
    file_biom = open(biom, "r")

    file_out1 = open(out_file1, "w")
    #file_out2 = open(out_file2, "w")

    line = file_biom.readline()
    dico_biom = json.loads(line.strip())
    print(dico_biom)
    # Pour la  1er collone : les cluster
    for i, value in enumerate(dico_biom["rows"]) :
        if dico_biom["rows"][i]["id"]:
            file_out1.write(dico_biom["rows"][i]["id"]+"\n")

    # Pour les entetes des autres collones 
    for i, value in enumerate(dico_biom["columns"]) :
        if dico_biom["columns"][i]["id"]:
            file_out1.write(dico_biom["columns"][i]["id"]+"\n")

    # Pour le contenu des collones (remplissage du tableau ) 
    for i, value in enumerate(dico_biom["data"]) :
        if dico_biom["data"][i][0]:
            print(["data"][i][1])
            #file_out1.write(dico_biom["data"][i][0]+"\n")

    file_biom.close()
    file_out1.close()

test(sys.argv[1], sys.argv[2])