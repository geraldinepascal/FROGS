import pandas as pd


def tsvParse():

    tsv_parse = "Reaction_prediction.tsv"
    #tsv_parse = "Markers"+"_"+args.output
    # je crée un dataframe
    df = pd.read_csv("/Users/moussa/FROGS_moussa/tools/FPStep2/DonneesPourTest/sorti_test_EC_etape2.tsv", header=0, delimiter='\t')
    # Je suprime la partie marker_count
    al = df.drop(columns=["metadata_NSTI", "16S_rRNA_Count"])
    print(df.drop(df.columns[[0, 1]], axis=1))
    print(al)
    # J écris dans un nouveau fichier
    al.to_csv(tsv_parse, "\t", index=None)

tsvParse()