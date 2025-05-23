#!/usr/bin/env Rscript

author = 'Ta Thi Ngan / Maria Bernard - Sigenae / Mahendra Mariadassou plateforme Migale'
copyright = 'Copyright (C) 2015 INRA'
license = 'GNU General Public License'
version = '1.1.0'
email = 'frogs-support@inrae.fr'
status = 'prod'

############## IMPORT

library(optparse)

############## OPTPARSE
option_list = list(
	make_option(c("-a", "--analysis"), type="character", default=NULL, help='Type of data to perform the differential analysis. ASV: DESeq2 is run on the ASVs abundances table. FUNCTION: DESeq2 is run on FROGSFUNC function abundances table (frogsfunc_functions_unstrat.tsv from FROGSFUNC function step).'),
	make_option(c("-i", "--inRdata"), type="character", default=NULL, help="The path of RData file containing a phyloseq object, result of FROGS Phyloseq Import Data (required)"),
	make_option(c("-v", "--var"), type="character", default=NULL, help="Experimental variable suspected to have an impact on ASVs abundances."),
	make_option(c("-f", "--inputFunction"), type="character", default=NULL, help='Input file of metagenome function predictions abundance. Required. (default: %(default)s).'),
	make_option(c("-s", "--samplefile"), type="character", default=NULL, help="path to sample file (format: TSV)."),
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
	version = paste( version, " [R : ",  R.version$major, ".",  R.version$minor, "; DESeq2 : ", as.character(DESeq2_version),"; Phyloseq : ", as.character(phyloseq_version), "]\n", sep="") 
	return(version)
}

if (opt$version){
	cat(get_version())
	quit()
}

########## MAIN
library(DESeq2)
library(phyloseq)

if (opt$analysis == "ASV"){
	load(opt$inRdata)
	cds <- phyloseq_to_deseq2(data, as.formula(paste("~", opt$var)))
	dds <- DESeq2::DESeq(cds, sfType = "poscounts")

}else if (opt$analysis == "FUNCTION"){
	sampleMetadata <- read.csv(file=opt$samplefile, sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
	# sampleMetadata <- sampleMetadata[ colnames(countData), , drop = FALSE]
	sampleMetadata[,opt$var] <- as.factor(sampleMetadata[,opt$var])
	
	inputFunction <- read.csv(file=opt$inputFunction, sep = '\t', header = TRUE, row.names = 3, check.names = FALSE)
	countData <- as.matrix(inputFunction[, c(4:ncol(inputFunction)) ])
	countData <- round(countData, 0)
	# select only samples with metadata
	countData <- countData[,rownames(sampleMetadata)]
	countData <- countData[!(rowSums(countData) == 0), !(colSums(countData) == 0)]
		# add log to check which samples / function have been removed.

	cds <- DESeq2::DESeqDataSetFromMatrix(countData, sampleMetadata, as.formula(paste("~",opt$var)))
	dds <- DESeq2::DESeq(cds)

}

save(dds, file=opt$outRdata)
