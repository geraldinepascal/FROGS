#!/bin/sh

#	to improve/simplify : il manque des noms pour les EC
#	le code pour KEGG n'a pas été développé pour FROGS et certaines parties ne sont pas nécessaire.
#	les noms level1 level2 level3 ne sont pas forcément appropriés pour COG (c'est plutôt categories, pathway, et pas de level3)

#	check si chaque colonne est unique ou redondante. est logique? Doit on imaginer une autre annotation?
#		les level sont utilisés pour créer un sunburst dans frogsfunc_functions. Est ce fonctionnel? Est ce pertinent ? Comme gérer la multiplicité des levels pour KEGG ou COG par exemple

# la concaténation des niveaux peut poser problème pour les stat


#################
# 							/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
#
# CONCLUSION : besoins de plus de réflexion ne communiquer que le fichier KEGG_hierarchy_all_ko_in_picrust à utiliser en dehors de FROGS
#
#
#################

#############################################################################################################################################
#												GET KEGG

#	some KO are not associated with pathway. NA fill each annotation level ==> check that everythin is ok in frogsfunc

# download kegg pathway name 
wget https://rest.kegg.jp/list/pathway
mv pathway KEGG_pathway.tsv

# download kegg hierarchy
wget https://rest.kegg.jp/get/br:ko00001
mv br\:ko00001 KEGG_hierarchy.txt

# download KEGG archive
wget https://www.genome.jp/ftp/db/kofam/archives/2019-03-20/ko_list.gz
gzip -d ko_list.gz 
mv ko_list KEGG_kofam_archive_2019_03_20.tsv

# parse KEGG_hierarchy and save Metabolism function in tsv file
python get_kegg_annot.py

# select picrust2 KEGG fun
gzip -dc ~/miniconda3/envs/frogsfunc@4.0.1/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz | wc -l
	# 13839 KEGG KO in picrust2

	# sur KO impliqué dans les metabolism uniquement
	gzip -dc ~/miniconda3/envs/frogsfunc@4.0.1/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz | awk -F "\t" 'BEGIN{
		while(getline<"KEGG_hierarchy_Metabolism_ko.tsv">0){tab[$1]=$0}
		}{
			if(tab[$1]!=""){
				print tab[$1]
			}
		}' > KEGG_hierarchy_Metabolism_ko_in_picrust.tsv
		
		wc -l KEGG_hierarchy_Metabolism_ko_in_picrust.tsv 
		# 13833
		grep -c "Metabolism;" KEGG_hierarchy_Metabolism_ko_in_picrust.tsv
		# 4218  
		grep -c "Removed Ortholog" KEGG_hierarchy_Metabolism_ko_in_picrust.tsv 
		# 519
		# 13833 -1 header - 519 removed = 13313 KO toujours dans KEGG
		# 13839 - 13832 = 7 KO picrust non retrouvé

	# sur all KO
	gzip -dc ~/miniconda3/envs/frogsfunc@4.0.1/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz | awk -F "\t" 'BEGIN{
		while(getline<"KEGG_hierarchy_all_ko.tsv">0){tab[$1]=$0}
		}{
			if(tab[$1]!=""){
				print tab[$1]
			}
		}' > KEGG_hierarchy_all_ko_in_picrust.tsv
		
		wc -l KEGG_hierarchy_all_ko_in_picrust.tsv 
		# 13833
		grep -c "Metabolism;" KEGG_hierarchy_all_ko_in_picrust.tsv
		# 4218  
		awk -F '\t' '$3==0 || $3=="NA"' KEGG_hierarchy_all_ko_in_picrust.tsv| wc -l 
		519  # 519 KO sans pathway
		grep -c "Removed Ortholog" KEGG_hierarchy_all_ko_in_picrust.tsv 
		# 519	==> en dehors de ces KO supprimés tous les autres sont associé à une hierarchy fonctionnelle
		# 13833 -1 header - 519 removed = 13313 KO 
		# 13839 - 13832 = 7 KO picrust non retrouvé


# check missing KO
gzip -dc ~/miniconda3/envs/frogsfunc@4.0.1/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/ko_info.tsv.gz | awk -F "\t" 'BEGIN{
	while(getline<"KEGG_hierarchy_Metabolism_ko.tsv">0){tab[$1]=$0}
	}{
		if(tab[$1]==""){print $0}
	}'

	K04225	TAKR; tachykinin-like receptor
	K06021	E3.6.3.27; phosphate-transporting ATPase [EC:3.6.3.27]
	K06022	E3.6.3.29; molybdate-transporting ATPase [EC:3.6.3.29]
	K06638	MAD1L; mitotic spindle assembly checkpoint protein MAD1
	K10826	E3.6.3.30; Fe3+-transporting ATPase [EC:3.6.3.30]
	K14049	DRD2N; dopamine D2-like receptor
	# set to probably removed KO ?

# format as gene_family_hierarchy.tsv
	# run KEGG_hierchy_reformat.Rmd


#############################################################################################################################################
#												GET COG

# download COG categories
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/fun-20.tab

# download COG description
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
mv cog-20.def.tab COG_descriptions.tsv
	# 1.	COG ID
	# 2.	COG functional category (could include multiple letters in the order of importance)
	# 3.	COG name
	# 4.	Gene associated with the COG (optional)
	# 5.	Functional pathway associated with the COG (optional)
	# 6.	PubMed ID, associated with the COG (multiple entries are semicolon-separated; optional)
	# 7.	PDB ID of the structure associated with the COG (multiple entries are semicolon-separated; optional)

# reformat as gene_family_hierarchy.tsv
	# run COG_annot.Rmd

#############################################################################################################################################
#												CREATE new gene_family_hierarchy.tsv

mv gene_family_hierarchy.tsv gene_family_hierarchy_frogs4.1.tsv
cp KEGG_annot.tsv gene_family_hierarchy.tsv 
tail -n +2 COG_annot.tsv >> gene_family_hierarchy.tsv 
awk -F "\t" '{if(substr($5,0,3) == "EC:"){print $0}}' gene_family_hierarchy_frogs4.1.tsv >> gene_family_hierarchy.tsv 

mv gene_family_hierarchy.tsv gene_family_hierarchy_KEGG_MetaboOnly.tsv
cp KEGG_all_annot.tsv gene_family_hierarchy.tsv 
tail -n +2 COG_annot.tsv >> gene_family_hierarchy.tsv 
awk -F "\t" '{if(substr($5,0,3) == "EC:"){print $0}}' gene_family_hierarchy_frogs4.1.tsv >> gene_family_hierarchy.tsv 