library(docopt)
library(future)
source("./functions_cite.R")
options(future.globals.maxSize = 64000 * 1024^4)
library("Seurat")
library("dplyr")
library("matrixStats")
library('tidyverse')
library(tidyseurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(sva)
library(harmony)
library(dsb)
library(docopt)
library(future)


doc <- "Usage:
  process_lanes.R [--rawdir=<rawdir>] [--demuxdir=<demuxdir>] [--save_name=<save_name>] [--lanes=<lanes>]
                  
Options:
  -h --help     Show this help message and exit.
  --lanes=<lanes>   Specify the range of lanes (e.g., 18:21).
  --rawdir=<rawdir>   Specify the output file name (e.g., output.RDS).
  --demuxdir=<demuxdir>
  --save_name=<save_name>
"

args <- docopt(doc)


#' Process Lanes Script
process_lanes <- function(lanes, rawdir, demuxdir, save_name) {
  
  
  timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
  
  # Define the name of the directory
  new_dir_name <- paste0(save_name,"_",timestamp)
  
  # Create the directory
  dir.create(new_dir_name)
  
  # Check if the directory was created successfully
  if (file.exists(new_dir_name) && file.info(new_dir_name)$isdir) {
    cat(paste("Directory '", new_dir_name, "' created successfully.\n", sep = ""))
  } else {
    cat(paste("Failed to create directory '", new_dir_name, "'.\n", sep = ""))
  }
  
  save_name <- paste0(save_name,"_",timestamp,"/",save_name,"_",timestamp)
  
B1_US_data <-read_h5_files(lanes, rawdir)
B1_US_SeuratObj <- create_seurat_objects2(B1_US_data, lanes)
##FOR CITE SEQ USE create_seurat_objects
#B1_US_SeuratObj <- create_seurat_objects(B1_US_data, lanes)
B1_US_merge <- merge_seurat_objects(B1_US_SeuratObj, lanes)
B1_US_merge <- add_metadata_and_filter(B1_US_merge, save_name)

}



process_lanes(as.numeric(unlist(strsplit(args$lanes, "[:,]"))), args$rawdir, args$demuxdir, args$save_name)

##Rscript ./process_downsampling.R --rawdir=./DATA/ --lanes=1,2,3,4,6,7,8 --save_name=SCALE


