#!/usr/bin/env Rscript

############## IMPORT

library(optparse)
library(DESeq2)
library(phyloseq)

############## MAIN
option_list = list(
	make_option(c("-i", "--inRdata"), type="character", default=NULL, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data (required)"),
	make_option(c("-v", "--var"), type="character", default=NULL, help="Experimental variable suspected to have an impact on OTUs abundances."),
	make_option(c("-o", "--outRdata"), type="character", default="DESeq2_preprocess.Rdata", help="The path to store resulting dataframe of DESeq2.[default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

opt

# check options
if (is.null(opt$inRdata)){
  print_help(opt_parser)
  stop("You need to provide one input file with the --in-Rdata\n", call.=FALSE)
}


load(opt$inRdata)

cds <- phyloseq_to_deseq2(data, as.formula(paste("~", opt$var)))

dds <- DESeq2::DESeq(cds, sfType = "poscounts")
save(dds, file=opt$outRdata)
