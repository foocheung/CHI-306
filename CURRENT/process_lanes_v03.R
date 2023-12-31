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
library(devtools)


# Define the command-line interface using docopt
doc <- "Usage:
  process_lanes.R [--rawdir=<rawdir>] [--demuxdir=<demuxdir>] [--save_name=<save_name>] [--lanes=<lanes>]
                  
Options:
  -h --help     Show this help message and exit.
  --lanes=<lanes>   Specify the range of lanes (e.g., 1,2,3,5).
  --rawdir=<rawdir>   Specify the input Directory (e.g., /DATA/).
  --demuxdir=<demuxdir> Specify the input DEMUX Directory (e.g., /DEMUX/).
  --save_name=<save_name> Specify the output Directory (e.g., OUTPUT).
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
  
  B1_US_data <- read_h5_files(lanes, rawdir, "23_306_3_")
  B1_US_SeuratObj <- create_seurat_objects(B1_US_data, lanes)
  B1_US_merge <- merge_seurat_objects(B1_US_SeuratObj, lanes)
  B1_US_merge <- add_metadata_and_filter(B1_US_merge, save_name)
  B1_US_merge <- read_demux_files(lanes, "S2_23_306_3_",demuxdir, B1_US_merge)
  data <- summarize_data(B1_US_merge)
  create_plots(data$lane_summ_dat, data$subj_summ_dat, save_name)
  B1_US_merge <- filter_cells_for_antibody_clumps(B1_US_merge)
  isotype.control.name.vec <- c(rownames(B1_US_merge[["CITE"]][c(61,63,66,117,115,97,137,83 ),]))
  B1_US_merge <- add_isotype_metadata(B1_US_merge,  isotype.control.name.vec)
  B1_US_merge <- filter_out_citeseq_outliers(B1_US_merge)
  #B1_US_merge <- DSB_Normalize_and_Save(B1_US_merge, isotype.control.name.vec)
  
  B1_US_merge <- normalize_and_scale_data(B1_US_merge)
  B1_US_merge <- run_rna_analysis(B1_US_merge)
  B1_US_merge <- run_cite_analysis(B1_US_merge)
  B1_US_merge <- run_wnn_analysis(B1_US_merge)
  B1_US_merge <- RunHarmony(B1_US_merge, group.by.vars = "BEST.GUESS")
  B1_US_merge <- run_rna_harmony(B1_US_merge)
  generate_umap_and_harmony_plots(B1_US_merge, save_name, lanes)
  B1_US_merge@meta.data$Lane<-as.factor(B1_US_merge@meta.data$Lane)
  save_seurat_as_rds(B1_US_merge, save_name)
}



process_lanes(as.numeric(unlist(strsplit(args$lanes, "[:,]"))), args$rawdir, args$demuxdir, args$save_name)
