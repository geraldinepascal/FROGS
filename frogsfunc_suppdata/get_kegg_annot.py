#!/usr/bin/env python3

# script download for an other project. Need to be updated to fit 


#####################################################################################################
# Téléchargement de l'arborescence kegg : https://rest.kegg.jp/get/br:ko00001
# Permet de récupérer les noms des fonctions
#	Fichier enregistré: KEGG_Orthology*.txt

#####################################################################################################
# to store KEGG, Brite, pathway class and class A and B names for class C
# and store function name

# Download pathway list and name https://rest.kegg.jp/list/pathway
# store which pathway are global map or overview (https://www.kegg.jp/kegg/pathway.html#metabolism  starint with 011 or 012)
d_pathAB_name_id = dict()
FH_pathAB_name_id = open("KEGG_pathway.tsv", "rt")

for line in FH_pathAB_name_id:
	(Id, Name) = line.strip().split("\t")
	suffix = Id.replace("map","")
	if suffix.startswith("011") or suffix.startswith("012"):
		Id = "ko"+suffix
		if suffix.startswith("011") :
			Class = "global map"
		else :
			Class = "overview map"
		d_pathAB_name_id[Id] = {"name" : Name,"class":Class}

d_brite_hierarchy=dict() # stocke les id, nom complet des pathway classA[;classB[;classC]]
d_fun_name_path = dict()	# stocke le nom et la liste d'ID de classe C

# parse pathway function composition hierarchy : https://rest.kegg.jp/get/br:ko00001
# format
	# #https://rest.kegg.jp/get/br:ko00001
	# +D	KO
	# !
	# A09100 Metabolism						# ==> équivaut à ko01100
	# B  09101 Carbohydrate metabolism
	# C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
	# D      K00844  HK; hexokinase [EC:2.7.1.1]
	# D      K12407  GCK; glucokinase [EC:2.7.1.2]
# FH_name = open("KEGG_Orthology_20231128.txt", "rt")
FH_name = open("KEGG_hierarchy.txt", "rt")
# is_in_metabolism = True
for line in FH_name:
	if not line[0] in ["A","B","C","D"]:
		continue
	if "B  09194 Poorly characterized" in line:		
		break
	# if "Genetic Information Processing" in line: 	# store only Metabolism function : 6232 KO sur / 20666 (ou 4218sur 13833 si on se retreint à celles dans picrust2)
	# 	is_in_metabolism = False

	if line[0] in ["A","B", "C"]:
		if line[0] == "A":
			line = line.replace("A","A ")
			line = line.replace("Brite Hierarchies","BRITE")
			line = line.replace("Not Included in Pathway or Brite","Other")

		(Class, Id, Name) = line.strip().split(maxsplit=2)

		if Class == "A" : 
			ClassA_name = Name

		if Class == "B" :
			ClassB_Id = Id
			ClassB_name = Name
			if ClassA_name == "BRITE":
				ClassB_name = ClassB_name.replace('Protein families: ','')
			if ClassA_name == "Other":
				ClassB_name = ClassB_name.replace('Unclassified: ','')

		if Class == "C" : 
			Name = Name.split("[")[0].strip()
			if "[PATH:ko" in line or "[BR:ko" in line:
				# ClassC_Id = "ko"+Id
				ClassC_Id = line.split()[-1].replace("[BR:","").replace("[PATH:","").replace("]","")
			else:
				ClassC_Id = Id

			if ClassB_Id  != "09112":
				# ignore 09112:
					# B  09112 Not included in regular maps
					# C    09113 Global maps only
					# D      K13997  PDHX; dihydrolipoamide dehydrogenase-binding protein of pyruvate dehydrogenase complex
				d_brite_hierarchy[ClassC_Id] = {"name" : ClassA_name + ";" + ClassB_name + ";" + Name, "class" : "C"}

	if line.startswith("D"):
		if ClassB_Id  == "09112" :
			# ignore 09112:
				# B  09112 Not included in regular maps
				# C    09113 Global maps only
				# D      K13997  PDHX; dihydrolipoamide dehydrogenase-binding protein of pyruvate dehydrogenase complex
			continue
		#D      K04363  PDGFRA, CD140A; platelet-derived growth factor receptor alpha [EC:2.7.10.1]
		d = line.strip().split(maxsplit=2)			#['D', 'K04363', 'PDGFRA, CD140A; platelet-derived growth factor receptor alpha [EC:2.7.10.1]']
		fun_id = d[1]

		if ";" in d[2]:
			d = d[2].split(";")						#['PDGFRA, CD140A','platelet-derived growth factor receptor alpha [EC:2.7.10.1]']
		else:
			d.pop(0)
		fun_name = d[1].strip()								#platelet-derived growth factor receptor alpha [EC:2.7.10.1]
		if not fun_id in d_fun_name_path:
			d_fun_name_path[fun_id] = {"name" : fun_name, "path" : list()}
		# if is_in_metabolism :
		d_fun_name_path[fun_id]["path"].append(ClassC_Id)
		# if fun_id == "K00046":
		# 	print(ClassA_name + ";" + ClassB_name + ";" + Name+ "\t"+str(is_in_metabolism))
		# if fun_id == "K00001":
		# 	print(ClassA_name + ";" + ClassB_name + ";" + Name)

# try to recover old function name : https://www.genome.jp/ftp/db/kofam/archives save in kofam_archives_2019-03-20_ko_list.tsv
FH_archives = open("KEGG_kofam_archive_2019_03_20.tsv")
for line in FH_archives:
	fun_id = line.split()[0]
	if fun_id != "knum" and not fun_id in d_fun_name_path:
		fun_name=line.split('\t')[-1].strip()
		d_fun_name_path[fun_id]= {"name" : fun_name, "path" : "Removed Ortholog"}


# print(d_fun_name_path["K00001"])

# faire un fichier bilan pour chaque fonction associer son nom le nombre de pathway KEGG, BRITE, OTHER (en 3 sections de colonnes séparées)
# FH_out = open("KEGG_hierarchy_Metabolism_ko.tsv" , "wt")			que Metabolism
FH_out = open("KEGG_hierarchy_all_ko.tsv" , "wt")
FH_out.write("Function_id\tFunction_name\tNumber_pathway\tPathway_ids\tPathway_names\n")
for fun in d_fun_name_path:
	if  d_fun_name_path[fun]["path"] == "Removed Ortholog":
		FH_out.write(fun + "\t" + d_fun_name_path[fun]["name"] + "\tNA\tNA\tRemoved Ortholog\n")
	else:
		path_id = list()
		path_name = list()
		for path in d_fun_name_path[fun]["path"]:
			if path in d_brite_hierarchy:
				path_id.append(path)
				path_name.append(d_brite_hierarchy[path]["name"])

		FH_out.write(fun + "\t" + d_fun_name_path[fun]["name"] + "\t" + str(len(d_fun_name_path[fun]["path"])) + "\t" + \
					  " & ".join(path_id) + "\t" + " & ".join(path_name) + "\n")

FH_out.close()