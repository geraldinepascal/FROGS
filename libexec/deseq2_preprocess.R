#!/usr/bin/env Rscript

############## IMPORT

library(optparse)

############## OPTPARSE
option_list = list(
	make_option(c("-i", "--inRdata"), type="character", default=NULL, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data (required)"),
	make_option(c("-v", "--var"), type="character", default=NULL, help="Experimental variable suspected to have an impact on OTUs abundances."),
	make_option(c("-o", "--outRdata"), type="character", default="DESeq2_preprocess.Rdata", help="The path to store resulting dataframe of DESeq2.[default= %default]"),
	make_option(c("--version"), action="store_true", default=FALSE, help="return version")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

############# VERSION
get_version <- function(){
	library(DESeq2)
	library(phyloseq)
	DESeq2_version = packageVersion("DESeq2")
	phyloseq_version = packageVersion("phyloseq")
	version = paste( " [R : ",  R.version$major, ".",  R.version$minor, "; DESeq2 : ", as.character(DESeq2_version),"; Phyloseq : ", as.character(phyloseq_version), "]\n", sep="") 
	return(version)
}

if (opt$version){
	cat(get_version())
	quit()
}

########## check options
if (is.null(opt$inRdata)){
  print_help(opt_parser)
  stop("You need to provide one input file with the --in-Rdata\n", call.=FALSE)
}

########## MAIN
library(DESeq2)
library(phyloseq)

load(opt$inRdata)

cds <- phyloseq_to_deseq2(data, as.formula(paste("~", opt$var)))

dds <- DESeq2::DESeq(cds, sfType = "poscounts")
save(dds, file=opt$outRdata)
