library(docopt)
library(future)
source("functions_cite.R")
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


Lanes <- c(1:16,18,19,21,22,23,24)
rawdir = "./306_2/"

files <- list()
for(i in 1:length(Lanes)){
  files[[i]] = paste(rawdir,"2multi_config_",
                     Lanes[[i]], "_sample_filtered_feature_bc_matrix.h5", sep="")
}

data <- list()
for(i in 1:length(Lanes)){
  data[[i]] = Read10X_h5(paste(rawdir,"2multi_config_",
                               Lanes[[i]],"_sample_filtered_feature_bc_matrix.h5", sep="")
  )
}

for(i in 1:length(Lanes)){
  rownames(data[[i]]$`Antibody Capture`) = gsub("\\.1$","", 
                                                rownames(data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
}


SeuratObj <- list()
for(i in 1:length(Lanes)){
  SeuratObj[[i]] <- CreateSeuratObject(counts =
                                         data[[i]]$'Gene Expression', assay = "RNA",
                                       min.feature = 20)
  
  
  if (i < 9){
  SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = data[[i]]$`Antibody Capture`[1:24,colnames(SeuratObj[[i]])])
  SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts =  data[[i]]$`Custom`[1:2,colnames(SeuratObj[[i]])])
  }
  else if (i > 8 & i < 22){
    SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = data[[i]]$`Antibody Capture`[1:27,colnames(SeuratObj[[i]])])
  }
    else{
      SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = data[[i]]$`Antibody Capture`[1:145,colnames(SeuratObj[[i]])])
      }
  
  
  
  SeuratObj[[i]] <- RenameCells(SeuratObj[[i]],
                                new.names = paste(substr(colnames(SeuratObj[[i]]),
                                                         start = 1, stop = 17),Lanes[[i]],"_306_2", sep = ""))
  SeuratObj[[i]]$Batch  <- rep("306_2",length(colnames(SeuratObj[[i]])))
  SeuratObj[[i]]$Lane <- rep(Lanes[i],length(colnames(SeuratObj[[i]])))
}


setwd("./DEMUX/")
demuxbestList = list()
for(i in 1:length(Lanes)){
  demuxbestList[[i]] = read.table(paste("S2_2multi_config_",
                                        Lanes[[i]],".best", sep = ""), sep = "\t", header = TRUE)
}
names(demuxbestList) <- Lanes

for(i in 1:length(Lanes)){
  demuxbestList[[i]]$NewBarcode = paste(substr(demuxbestList[[i]]$BARCODE, 
                                               start = 1, stop = 17),
                                        Lanes[[i]],"_306_2", sep = "")
}

#check the overlay of cells between seuratobj and demux list
for(i in 1:length(SeuratObj)){
  length(which(colnames(SeuratObj[[i]]) %in% demuxbestList[[i]]$NewBarcode))
  print(setdiff(colnames(SeuratObj[[i]]), demuxbestList[[i]]$NewBarcode))
}


demuxbestdf <- plyr::ldply(demuxbestList, data.frame)
rownames(demuxbestdf) <- demuxbestdf$NewBarcode

for(i in 1:length(SeuratObj)){
  SeuratObj[[i]]  <- AddMetaData(SeuratObj[[i]],
                                 metadata = demuxbestdf[colnames(SeuratObj[[i]]),])
}


ObjList<-SeuratObj


ObjList_neg = list()
ObjList_SNG = list()
for(i in 1:length(ObjList)){
  ObjList_neg[[i]] <- subset(ObjList[[i]], DROPLET.TYPE == "AMB"|nFeature_RNA < 100)
  # ObjList_neg[[i]] <- subset(ObjList[[i]], cells = sample(Cells(ObjList_neg[[i]]),1000))
  ObjList_SNG[[i]] <- subset(ObjList[[i]], nFeature_RNA > 100 & DROPLET.TYPE == "SNG")
}



for(i in 1:length(ObjList_SNG)){
  ObjList_SNG[[i]][["percent.mito"]] <- PercentageFeatureSet(object = ObjList_SNG[[i]], 
                                                             pattern = "^MT-")
  }


for(i in 1:length(ObjList_SNG)){
  ObjList_SNG[[i]] = subset(ObjList_SNG[[i]], subset = nFeature_RNA > 200 &nFeature_RNA < 5000 & nCount_RNA < 30000 & percent.mito < 30 & nCount_CITE < 20000)
}



VlnPlot(ObjList_SNG[[1]], features= c("nCount_RNA","nCount_CITE", "nFeature_RNA","percent.mito"), 
        group.by = "Lane",pt.size = 0, ncol = 2)

merged_obj <- merge(ObjList_SNG[[1]],ObjList_SNG[2:length(ObjList_SNG)])
VlnPlot(merged_obj, features= c("nCount_RNA","nCount_CITE", "nFeature_RNA","percent.mito"), 
        group.by = "Lane",pt.size = 0, ncol = 2)
ggsave("merged.pdf")

saveRDS(merged_obj,"306_2_merge_with_SNP_v1.rds")

saveRDS(ObjList_SNG,"306_2_ObjList_SNG_with_SNP_v1.rds")
