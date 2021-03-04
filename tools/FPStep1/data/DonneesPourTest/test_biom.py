import os
import sys
import json
import re
def test(file_tree, biom, out_file1, out_file2):
    file = open(file_tree, "r")
    line = file.readline()
    #List of cluster
    list_cluster = []
    #Boucle sur le fichier tree
    #Je splite sur (,) , regex sur le cluster et récupérer le groupe1
    #Je parcour ligne par ligne
    while line: 
        for i, v in enumerate(line.split(",")):
            group = re.search("(Cluster_[0-9]+)", v)
            if group:
                ide = group.group(1)
                list_cluster.append(ide)

        line = file.readline() 

    file.close()

    file_biom = open(biom, "r")

    file_out1 = open(out_file1, "w")
    file_out2 = open(out_file2, "w")

    line = file_biom.readline()
    dico_biom = json.loads(line.strip())
    print(dico_biom)

    for i, value in enumerate(dico_biom["rows"]) :
        if dico_biom["rows"][i]["id"] in list_cluster :
            file_out1.write(dico_biom["rows"][i]["id"]+"\n")
        else:
            file_out2.write(dico_biom["rows"][i]["id"]+"\n")


    file_biom.close()
    file_out1.close()
    file_out2.close()



test(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

