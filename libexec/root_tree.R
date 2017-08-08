#!/usr/bin/env Rscript
##################################

author = 'Maria Bernard - Sigenae'
copyright = 'Copyright (C) 2015 INRA'
license = 'GNU General Public License'
version = '1.0.0'
email = 'frogs@inra.fr'
status = 'prod'


############ IMPORT ###############
library(ape)
library(phangorn)

########### FUNCTIONS #############

usage <- function(){
# 
message <- paste(sep="\n" ,"# Description : ", "\troot_tree.R is dedicated to root tree thanks to a midpoint method",
	paste("\troot_tree.R",version ),
	"# Usage : ","\tRscript root_tree.R input_tree.nwk output_rooted_tree.nwk",
	"",
	"  input_tree.nwk : \t\t Path file of input tree to root" ,
	"  output_rooted_tree.nwk : \t Path file of output rooted tree to save",
	"")
	
return(message)
}

root_tree <- function(in_tree, out_tree){
	tree <- read.tree(in_tree)
	rooted_tree <- phangorn::midpoint(tree)
	write.tree(rooted_tree, file=out_tree)
}


########### MAIN #################
# Get all arguments
args <- commandArgs(trailingOnly = TRUE)

# check params
if (args[1]=="-h" || length(args)==0){
	cat(usage())
	quit()
}
if (args[1]=="-v"){
	cat("root_tree.R : ", version ,"\n") 
	quit()
}
if(length(args) != 2 ){
 	stop(paste("You provide ", length(args), " arguments instead of 2 asked!", "\n\n",usage()))
}

# process 
cat("Input tree file:", args[1], "\n")
cat("Output rooted tree file:", args[2], "\n")

root_tree(args[1], args[2])
