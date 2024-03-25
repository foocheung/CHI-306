library(docopt)
library(future)
#source("/Users/xuq6/Documents/Vaccination project/CITEseq data/306-3/functions_cite.R")
source("306_4_functions_cite.R")
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
  process_lanes.R [--rawdir=<rawdir>] [--demuxdir=/Volumes/LISB_CSI/306_MERGED_RUN1_RUN2/DEMUX] [--save_name=/Users/xuq6/Documents/Vaccination project/CITEseq data/306-4/20210116 Out] [--lanes=1:24]
                  
Options:
  -h --help     Show this help message and exit.
  --lanes=<lanes>   Specify the range of lanes (e.g., 1,2,3,5).
  --rawdir=<rawdir>   Specify the input Directory (e.g., /DATA/).
  --demuxdir=<demuxdir> Specify the input DEMUX Directory (e.g., /DEMUX/).
  --save_name=<save_name> Specify the output Directory (e.g., OUTPUT).
"

args <- docopt(doc)


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
  
  B1_US_data <- read_h5_files(lanes, rawdir, "multi_config_")
  
  B1_US_SeuratObj <- create_seurat_objects(B1_US_data, lanes)
#  B1_US_merge <- merge_seurat_objects(B1_US_SeuratObj, lanes)
  B1_US_merge <- read_demux_files(lanes, "S2_multi_config_",demuxdir, B1_US_SeuratObj)
  
  SeuratObj<-B1_US_merge
  
  # library(dsb)
   Spos_SeuratObj <-list(SeuratObj[[4]],SeuratObj[[12]],SeuratObj[[20]])
   Sneg_SeuratObj <- merge(SeuratObj[[1]],c(SeuratObj[2:3],SeuratObj[9:11],SeuratObj[17:19]))
   CD4_SeuratObj <-merge(SeuratObj[[5]],c(SeuratObj[6],SeuratObj[13:14],SeuratObj[21:22]))
   CD8_SeuratObj <-merge(SeuratObj[[7]],c(SeuratObj[8],SeuratObj[15:16],SeuratObj[23:24]))
   B_SeuratObj <- merge(Sneg_SeuratObj,Spos_SeuratObj)
   
   ObjList <- list(B_SeuratObj,CD4_SeuratObj,CD8_SeuratObj)

  rm(list=c("Spos_SeuratObj","Sneg_SeuratObj","CD4_SeuratObj","CD8_SeuratObj","B_SeuratObj","B1_US_merge"))
  ObjList<- filter_cells_for_antibody_clumps(ObjList)
  
  
##  saveRDS(ObjList,paste(new_dir_name,"/","ObjList.rds", sep=""))
  
  ObjList_filt<- DSB_Normalize(ObjList)
  
  rm("ObjList")
 


saveRDS(ObjList_filt,paste(new_dir_name,"/","ObjList_filt.rds", sep=""))
 
  names(ObjList_filt) <- c("Bcells","CD4","CD8")
  
  ObjList_filt<- filter(ObjList_filt)
  merged_obj <- merge(ObjList_filt[[1]],ObjList_filt[2:3])
  VlnPlot(merged_obj, features= c("nCount_RNA","nCount_CITE", "nFeature_RNA","percent.mito"), 
          group.by = "Lane",pt.size = 0, ncol = 2)
  ggsave(paste(new_dir_name,"/","QC.pdf",sep=""),width =15, height = 10)
  Idents(ObjList_filt[[1]]) <- "Lane"
  
  metadata <- merged_obj@meta.data
  colnames(metadata)
  table(metadata$DROPLET.TYPE)


##make function here

  summary_data <- metadata %>% group_by(Lane)%>%
    summarise(mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
              mean_nCount_CITE = mean(nCount_CITE, na.rm = TRUE),
              mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
              mean_percent.mito= mean(percent.mito, na.rm = TRUE))
  
  write.csv(summary_data,paste(new_dir_name,"/","summary_data.csv",sep=""))
 

  sng_data <-metadata %>%
  group_by(Lane,DROPLET.TYPE) %>%
  summarise(n = n())
  
  write.csv(sng_data,paste(new_dir_name,"/","sng_data.csv",sep=""))
  
  
 
  
  ##Rscript ./306_4_process.R --rawdir=./CHI_4_Merged/ --demuxdir=./DEMUX/ --save_name=TEST  --lanes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
  
  
  
  
    
}

process_lanes(as.numeric(unlist(strsplit(args$lanes, "[:,]"))), args$rawdir, 
              args$demuxdir, args$save_name)
  
  



