import os
import sys
import json
import re
def excluded_sequence(file_tree, file_fasta, out_file):

    """
    @summary: Returns the excluded sequence.
    @param fasta_file: [str] Path to the fasta file to process.
    @param tree_file: [str] Path to the tree file to process.
    @return: [int] The file of no aligned sequence.
    """
    #Lecture du fichier tree
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

    file_fasta = open(file_fasta, "r")
    file_out   = open(out_file, "w")

    line = file_fasta.readline()
    print(list_cluster)

    while line:
        if line[0] == ">":
            if line[1:].strip() not in list_cluster:
                if f < 2:
                    print(line.strip())
                file_out.write(line.strip()+"\n")
                line = file_fasta.readline()
                file_out.write(line.strip()+"\n")


        line = file_fasta.readline()

    file_fasta.close()
    file_out.close()

excluded_sequence(sys.argv[1], sys.argv[2], sys.argv[3])
