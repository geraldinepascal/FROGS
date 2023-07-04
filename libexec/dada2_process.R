#!/usr/bin/env /usr/local/public/R/bin/Rscript

author = 'Olivier Rué'
copyright = 'Copyright (C) 2023 INRAE'
license = 'GNU General Public License'
version = '1.0'
email = 'frogs-support@inrae.fr'
status = 'dev'

############## IMPORT

library(optparse)

############## OPTPARSE
option_list = list(
  make_option(c("--R1Files"), type="list", default=NULL, help="List of R1 files to be process"),
  make_option(c("--R2Files"), type="list", default=NULL, help="List of R2 files to be process"),
  make_option(c("--pseudopooling"), action="store_true", default=FALSE, help="Perform pseudo-pooling to reduce inconvenients of independent sample processing"),
  make_option(c("-o", "--outputDir"), type="character", default=".", help="The directory path to write denoised FASTQ files. [default= %default]"),
  make_option(c("-f", "--fileNames"), type="character", default=".", help="Linked R1 and R2 files in the current analysis."),
  make_option(c("-t", "--threads"), type="integer", default=1, help="Number of CPUs to use. [default= %default]"),
  make_option(c("--version"), action="store_true", default=FALSE, help="return version"),
  make_option(c("--debug"), action="store_true", default=FALSE, help="debug mode")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

############# VERSION
get_version <- function(){
  library(dada2)
  dada2_version = packageVersion("dada2")
  version = paste( version, " [R : ",  R.version$major, ".",  R.version$minor, "; dada2 : ", as.character(dada2_version),"]\n", sep="") 
  return(version)
}

if (opt$version){
  cat(get_version())
  quit()
}

########## check options
if (is.null(opt$R1Files)){
  print_help(opt_parser)
  stop("You need to provide the list of R1 FASTQ files to be denoised with the --R1Files parameter\n", call.=FALSE)
}

########## Functions
getDerep <- function(object, ...) {
  if(is(object, "derep")) { object }
  else if(is.list.of(object, "derep")) { object }
  else if(is(object, "character")) { derepFastq(object, ...) }
  else{
    stop("Unrecognized format: Requires derep-class object, list of derep-class objects, a directory containing fastq files, or a character vector of fastq files.")
  }
}

writeFASTQsingle <- function(data, file, direction = c("forward", "reverse")) {
  direction <- match.arg(direction)
  string <- paste0("@", data$id, ";size=", data$abundance, "\n", data[[direction]], "\n+\n", 
                   data[[paste0(direction, "Qual")]])
  writeLines(string, con = file)
}

is.list.of <- function(x, ctype) {
  if(!is.list(x)) return(FALSE)
  else return(all(sapply(x, is, ctype)))
}

write2FastqFromDada <- function(dadaF, derepF, dadaR, derepR, path)
{
  nrecs <- c(length(dadaF), length(derepF), length(dadaR), 
             length(derepR))
  if (length(unique(nrecs)) > 1) 
    stop("The dadaF/derepF/dadaR/derepR arguments must be the same length.")
  rval <- lapply(seq_along(dadaF), function(i) {
    mapF <- getDerep(derepF[[i]])$map
    mapR <- getDerep(derepR[[i]])$map
    if (!(is.integer(mapF) && is.integer(mapR))) 
      stop("Incorrect format of $map in derep-class arguments.")
    if (!(length(mapF) == length(mapR) && max(mapF, na.rm = TRUE) == 
          length(dadaF[[i]]$map) && max(mapR, na.rm = TRUE) == 
          length(dadaR[[i]]$map))) {
      stop("Non-corresponding derep-class and dada-class objects.")
    }
    rF <- dadaF[[i]]$map[mapF]
    rR <- dadaR[[i]]$map[mapR]
    pairdf <- data.frame(id = "", abundance = 0, forward = rF, 
                         reverse = rR)
    ups <- unique(pairdf)
    keep <- !is.na(ups$forward) & !is.na(ups$reverse)
    ups <- ups[keep, ]
    if (nrow(ups) == 0) {
      outnames <- c("abundance", "forward", 
                    "reverse")
      ups <- data.frame(matrix(ncol = length(outnames), 
                               nrow = 0))
      names(ups) <- outnames
      
      return(ups)
    }
    else {
      int_to_quality <- function(qualities) {
        qualities <- qualities[!is.na(qualities)]
        intToUtf8(as.integer(floor(qualities)) + 33L)
      }
      Funqseq <- unname(as.character(dadaF[[i]]$clustering$sequence[ups$forward]))
      Runqseq <- unname(as.character(dadaR[[i]]$clustering$sequence[ups$reverse]))
      Funqqual <- apply(dadaF[[i]]$quality[ups$forward, , drop = FALSE], MARGIN = 1, int_to_quality)
      Runqqual <- apply(dadaR[[i]]$quality[ups$reverse, , drop = FALSE], MARGIN = 1, int_to_quality)
      
      tab <- table(pairdf$forward, pairdf$reverse)
      ups$abundance <- tab[cbind(ups$forward, ups$reverse)]
      
      #ups <- ups[order(ups$abundance, decreasing = TRUE), ]
      
      rownames(ups) <- NULL
      ups$id <- 1:nrow(ups)
      ups$forward <- Funqseq
      ups$reverse <- Runqseq
      ups$forwardQual <- Funqqual
      ups$reverseQual <- Runqqual
      sample_name <- names(dadaF)[i]
      R1_path <- file.path(path,paste0(sample_name,"_denoised_R1.fastq"))
      R2_path <- file.path(path,paste0(sample_name,"_denoised_R2.fastq"))
      writeFASTQsingle(ups, file=R1_path, direction="forward")
      writeFASTQsingle(ups, file=R2_path, direction="reverse")
      #set up writing
      cat(R1_path, R2_path, file=opt$fileNames, append=TRUE, sep = ",")
      cat("", file=opt$fileNames, append=TRUE, sep = "\n")

      return(ups)
    }
  })
  if (!is.null(names(dadaF))) 
    names(rval) <- names(dadaF)
  if (length(rval) == 1) 
    rval <- rval[[1]]
  return(rval)
}

write1FastqFromDada <- function(dadaF, derepF, path)
{
  if (is(dadaF, "dada")) 
    dadaF <- list(dadaF)
  if (is(derepF, "derep")) 
    derepF <- list(derepF)
  else if (is(derepF, "character") && length(derepF) == 1 && 
           dir.exists(derepF)) 
    derepF <- parseFastqDirectory(derepF)
  if (!(is.list.of(dadaF, "dada"))) {
    stop("dadaF must be provided as dada-class objects or lists of dada-class objects.")
  }
  if (!((is.list.of(derepF, "derep") || is(derepF, "character")))) {
    stop("derepF must be provided as derep-class objects or as character vectors of filenames.")
  }
  nrecs <- c(length(dadaF), length(derepF))
  if (length(unique(nrecs)) > 1) 
    stop("The dadaF/derepF arguments must be the same length.")
    
    
  rval <- lapply(seq_along(dadaF), function(i) {
    mapF <- getDerep(derepF[[i]])$map
    if (!(is.integer(mapF))) 
      stop("Incorrect format of $map in derep-class arguments.")
    if (!(max(mapF, na.rm = TRUE) == 
          length(dadaF[[i]]$map) )) {
      stop("Non-corresponding derep-class and dada-class objects.")
    }
    rF <- dadaF[[i]]$map[mapF]
    pairdf <- data.frame(id = "", abundance = 0, forward = rF)
    ups <- unique(pairdf)
    keep <- !is.na(ups$forward)
    ups <- ups[keep, ]
    if (nrow(ups) == 0) {
      outnames <- c("abundance", "forward")
      ups <- data.frame(matrix(ncol = length(outnames), 
                               nrow = 0))
      names(ups) <- outnames
      if (verbose) {
        message("No paired-reads (in ZERO unique pairings) successfully merged out of ", 
                nrow(pairdf), " pairings) input.")
      }
      return(ups)
    }
    else {
      int_to_quality <- function(qualities) {
        qualities <- qualities[!is.na(qualities)]
        intToUtf8(as.integer(floor(qualities)) + 33L)
      }
      Funqseq <- unname(as.character(dadaF[[i]]$clustering$sequence[ups$forward]))
      Funqqual <- apply(dadaF[[i]]$quality[ups$forward, , drop = FALSE], MARGIN = 1, int_to_quality)
      
      
      tab <- table(pairdf$forward)
      ups$abundance <- tab[cbind(ups$forward)]
      
      ups <- ups[order(ups$abundance, decreasing = TRUE), 
      ]
      
      rownames(ups) <- NULL
      ups$id <- 1:nrow(ups)
      ups$forward <- Funqseq
      ups$forwardQual <- Funqqual
      sample_name <- names(dadaF)[i]
      R1_path <- file.path(path,paste0(sample_name,"_denoised_R1.fastq"))
      writeFASTQsingle(ups, file=R1_path, direction="forward")
      #set up writing
      cat(R1_path, file=opt$fileNames, append=TRUE, sep = ",")
      cat("", file=opt$fileNames, append=TRUE, sep = "\n")

      return(ups)
    }
  })
  if (!is.null(names(dadaF))) 
    names(rval) <- names(dadaF)
  if (length(rval) == 1) 
    rval <- rval[[1]]
  return(rval)
}

########## MAIN
if(file.exists("dadaFs.rds") && file.exists("dadaRs.rds") && file.exists("derepFs.rds") && file.exists("derepRs.rds")){
	dadaFs <- readRDS("dadaFs.rds")
	dadaRs <- readRDS("dadaRs.rds")
	derepFs <- readRDS("derepFs.rds")
	derepRs <- readRDS("derepRs.rds")
}else{

library(dada2)

# store R1 files
fnFs <- strsplit(opt$R1Files, ",")[[1]]

# store R2 files
saveRDS(fnFs,"fnFs.rds")
if(!is.null(opt$R2Files)){
	fnRs <- strsplit(opt$R2Files, ",")[[1]]
}

saveRDS(fnRs,"fnRs.rds")
# function to get samples from file names
get.sample.name <- function(fname) paste(strsplit(basename(fname), "_R1.fastq.gz")[[1]][1],collapse="_")
# get sample names
sample.names <- unname(sapply(fnFs, get.sample.name))
saveRDS(sample.names,"samples.rds")

### Learn the Error Rates
errF <- learnErrors(fnFs, multithread=opt$threads)
if(!is.null(opt$R2Files)){
	errR <- learnErrors(fnRs, multithread=opt$threads)
}
###

### Dereplicate 
derepFs <- derepFastq(fnFs, verbose = TRUE)
names(derepFs) <- sample.names
if(!is.null(opt$R2Files)){
	derepRs <- derepFastq(fnRs, verbose = TRUE)
	names(derepRs) <- sample.names
}

### Sample Inference
if(opt$pseudopooling){
	dadaFs <- dada(derepFs, err=errF, multithread=opt$threads, pool=TRUE)
	if(!is.null(opt$R2Files)){
		dadaRs <- dada(derepRs, err=errF, multithread=opt$threads, pool=TRUE)
	}
}else{
	dadaFs <- dada(derepFs, err=errF, multithread=opt$threads)
	if(!is.null(opt$R2Files)){
		dadaRs <- dada(derepRs, err=errF, multithread=opt$threads)
	}
}

if (opt$debug){
    saveRDS(derepFs, file.path(opt$outputDir,"derepFs.rds"))
    saveRDS(dadaFs, file.path(opt$outputDir,"dadaFs.rds"))
    if(!is.null(opt$R2Files)){
		saveRDS(derepRs, file.path(opt$outputDir,"derepRs.rds"))
		saveRDS(dadaRs, file.path(opt$outputDir,"dadaRs.rds"))
	}
}

if(!is.null(opt$R2Files)){
	write2FastqFromDada(dadaFs, derepFs, dadaRs, derepRs, path=opt$outputDir)
}else{
	write1FastqFromDada(dadaFs, derepFs, path=opt$outputDir)
}
}
