#!/usr/bin/env python3
import pandas as pd

#df.drop(['B', 'E'], axis='columns', inplace=True)
# or df = df.drop(['B', 'E'], axis=1) without the option inplace=True

tsv_parse = "resultat_parse_KO.tsv"
#tsv_parse = "DonneesPourTest/sorti_etape2.tsv"
    # je crée un dataframe
#df = pd.read_csv("DonneesPourTest/sorti_etape2.tsv", header=0, delimiter='\t')
# Pour EC
#df = pd.read_csv("../hsp/DonneesPourTest/sorti_test1_etape2.tsv", header=0, delimiter='\t')
# Pour KO
df = pd.read_csv("../hsp/DonneesPourTest/sorti_test1_KO_etape2.tsv", header=0, delimiter='\t')

# Pour Tout
#df = pd.read_csv("out_test.tsv", header=0, delimiter='\t')

    # Je suprime les mauvaises colonnes
#al = df.drop(['metadata_NSTI', '16S_rRNA_Count'], axis='columns', inplace=True)

al = df.drop(columns=["metadata_NSTI", "16S_rRNA_Count"])
#df.drop(columns=["B", "C"])
print(al)

al.to_csv(tsv_parse, "\t", index=None)

print(tsv_parse)


#################################





#df.drop(['B', 'E'], axis='columns', inplace=True)
    #al = df.drop(df.columns[[0, 1, 3]], axis=1)
    #print(df.drop(df.columns[[0, 1]], axis=1))
    #print(al)

# def tsvParse(out_tsv):

#     tsv_parse = "out_tsv.tsv"
#     # je crée un dataframe
#     df = pd.read_csv(out_tsv, header=0, delimiter='\t')
#     # Je suprime les mauvaises colonnes
#     al = df.drop(df.columns[[0, 1, 3]], axis=1)
#     #print(df.drop(df.columns[[0, 1]], axis=1))
#     print(al)
#     # J écris dans un nouveau fichier
#     al.to_csv(tsv_parse, "\t", index=None)


#     metadata_NSTI	16S_rRNA_Count