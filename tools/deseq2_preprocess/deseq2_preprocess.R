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
data_deseq2 <- data

## Add 1 to all count to eliminite huge number of 0
otu_table(data_deseq2) <- otu_table(data) + 1

import_cmd <- paste('cds <- phyloseq_to_deseq2(data_deseq2, ~', opt$var, ')')
eval(parse(text = import_cmd))

dds <- DESeq2::DESeq(cds)
save(dds, file=opt$outRdata)
