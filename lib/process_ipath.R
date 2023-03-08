#!/usr/bin/env Rscript
##################################

author = 'Vincent Darbot - GenPhySe'
copyright = 'Copyright (C) 2023 INRA'
license = 'GNU General Public License'
version = '1.0.0'
email = 'frogs-support@inrae.fr'
status = 'prod'

############ IMPORT ###############
library(optparse)
library(dplyr)
library(httr)
library(stringr)
########### ARGS #############

########### FUNCTIONS #############

generate_ipath_input <- function(file) {
  out_str = ""
  FH_in <- readLines(file)
  for (li in FH_in) {
    if (startsWith(li, "FUNC")) {
      next
    } else {
      li <- unlist(strsplit(li, " "))
      out_str <- paste0(out_str, paste(li[1], li[2], li[3], sep = " "), "\n")
    }
  }
  return(out_str)
}

to_parameters <- function(
    ipath_data ,
    export_type='svg',
    include_metabolic=TRUE,
    include_secondary=FALSE,
    include_antibiotic=FALSE,
    include_microbial=FALSE,
    whole_modules=FALSE,
    whole_pathways=FALSE,
    keep_colors=FALSE,
    default_opacity=1,
    default_width=3,
    default_radius=7,
    default_color='#c7c5c5',
    query_reactions=FALSE,
    tax_filter='',
    export_dpi=1200) {

      # only svg working
      # allowed_export_types= list('svg','png','pdf','eps')

      ipath_parameters =  list(selection=ipath_data,
        export_type=export_type,
        keep_colors=ifelse(keep_colors, 1, 0),
        include_metabolic=ifelse(include_metabolic, 1, 0),
        include_secondary= ifelse(include_secondary, 1, 0),
        include_antibiotic= ifelse(include_antibiotic, 1, 0),
        include_microbial= ifelse(include_microbial, 1, 0),
        whole_modules= ifelse(whole_modules, 1, 0),
        whole_pathways=ifelse(whole_pathways, 1, 0),
        default_opacity= default_opacity,
        default_width= default_width,
        default_color= default_color,
        default_radius= default_radius,
        query_reactions= query_reactions,
        tax_filter= tax_filter,
        export_dpi=export_dpi
      )
      return(ipath_parameters)
}

get_map <- function(selection, map_name = "map") {
    library(httr)
    url <- "https://pathways.embl.de/mapping.cgi"

    parameters <- to_parameters(selection)
    r <- POST(url, body = parameters)
    stop_for_status(r)
    text <- content(r, "text")
    write(text, file = paste0(map_name))
}


## MAIN
process_iPath <- function(input_file, output_file) {
    func_str <- generate_ipath_input(input_file)
    get_map(func_str, output_file)
}
