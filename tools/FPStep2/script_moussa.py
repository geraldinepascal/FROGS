#!/usr/bin/env python3
import json
from collections import OrderedDict

dict_tsv = OrderedDict()

FH = open("fichier_0_out_test.tsv","rt")

keys = FH.readline().strip().split()

for line in FH:   
    for idx,value in enumerate(line.split()):
        if idx == 0:
            seq_id = value
            if seq_id not in dict_tsv:
                dict_tsv[seq_id] = OrderedDict()
        else :
            annot = keys[idx]
            dict_tsv[seq_id][annot] = value

FH.close()

# FH = open("fichier_1_out_test.tsv","rt")

# keys = FH.readline().strip().split()

# for line in FH:   
#     for idx,value in enumerate(line.split()):
#         if idx == 0:
#             seq_id = value
#             if seq_id not in dict_tsv:
#                 dict_tsv[seq_id] = OrderedDict()
#         else :
#             annot = keys[idx]
#             dict_tsv[seq_id][annot] = value            

print(dict_tsv)


# fout = open("Mousse.txt", "w")
# print("test1")
# # for k, v in dict_tsv.items():
# #     fout.write(str(k) + "\t" + str(v) + "\n" )
# # fout.close()
# header ="sequence"
# for seq_id in dict_tsv:
#     if header == "sequence":
#         header = header + "\t" + "\t".join(dict_tsv[seq_id].keys()) + "\n"
#         fout.write(header)
#     line=seq_id
#     for annot in dict_tsv[seq_id]:
#         line=line + "\t" + dict_tsv[seq_id][annot]
#     fout.write(line + "\n")  
# fout.close()
# ##############
# print("test2")

