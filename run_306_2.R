# Load necessary libraries
library("Seurat")
library("dplyr")
library("matrixStats")
library('tidyverse')
library(tidyseurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(docopt)
library(harmony)
library(SingleR)

source("./functions_cite.R")

#Rscript script.R --lanes=18,19,20,21 --output=v2306_2_seur_18_21.RDS
# Define the command line arguments using docopt
doc <- "
Usage:
  script.R [--lanes=<lanes>] [--output=<output>]

Options:
  -h --help     Show this help message and exit.
  --lanes=<lanes>   Specify the range of lanes (e.g., 18:21).
  --output=<output>   Specify the output file name (e.g., output.RDS).
"

# Parse the command line arguments
args <- docopt(doc)

# Extract arguments
lanes <- args$lanes
output <- args$output

# Convert lanes to a numeric vector
B1LanesUnsorted <- as.numeric(unlist(strsplit(lanes, "[:,]")))


B1_US_data = list()


rawdir = "./Ver2_RUN_306_2/DATA/"


for(i in 1:length(B1LanesUnsorted)){
  B1_US_data[[i]] = Read10X_h5(paste(rawdir,"23_306_2_",B1LanesUnsorted[i],"_sample_filtered_feature_bc_matrix.h5", sep=""))
}


B1_US_SeuratObj = list()


for(i in 1:length(B1LanesUnsorted)){
  B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = (B1_US_data[[i]]$'Gene Expression'), assay = "RNA", min.feature = 0)
  B1_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:142,])
  B1_US_SeuratObj[[i]][["Custom"]] <- CreateAssayObject(counts=B1_US_data[[i]]$'Custom'[1:3,])
  B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17),B1LanesUnsorted[[i]], sep = ""))
  B1_US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_US_SeuratObj[[i]])))
  B1_US_SeuratObj[[i]]$Lane  <- rep(B1LanesUnsorted[[i]], length(colnames(B1_US_SeuratObj[[i]])))
  
}


names(B1_US_data) = names(B1_US_SeuratObj) = B1LanesUnsorted

B1_US_merge = merge(B1_US_SeuratObj[[1]], B1_US_SeuratObj[2:length(B1LanesUnsorted)])
mito.genes = grep(pattern = "^MT-", x = rownames(B1_US_merge), value = TRUE)
B1_US_merge<-AddMetaData(object = B1_US_merge, metadata = Matrix::colSums(B1_US_merge[mito.genes,])/Matrix::colSums(B1_US_merge), col.name = "percent.mito")


# df<-B1_US_merge %>%
#   group_by(Lane) %>%
#   summarise(across(nCount_RNA, median), cell_count = n()) %>%
#   mutate("Comment"="Run1: Filtered, same results as cell ranger")
# 
# 
# df2<-B1_US_merge %>%
#   group_by(Lane) %>% filter(nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 20000 & percent.mito < 0.30) %>%
#   summarise(across(nCount_RNA, median), cell_count = n()) %>%
#   mutate("Comment"="RUN2: Same as Run1 + nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 20000 & percent.mito < 0.30")
# 
# 
# 
# all<-rbind(df,df2)
# 
# all<-janitor::clean_names(all) %>% as.matrix()
# 
# write.table(all %>% as.matrix(),"all.txt")


pdf("prefiltering.pdf", height=21)
VlnPlot(B1_US_merge, c( "percent.mito","percent_ribo","percent_hb", "nFeature_RNA", "nCount_RNA"), group.by = "Lane", pt.size=0, ncol = 1)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

B1_US_merge<-B1_US_merge %>%
  filter(nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 20000 & percent.mito < 0.30)

B1_US_merge <- PercentageFeatureSet(B1_US_merge, "^RP[SL]", col.name = "percent_ribo")
B1_US_merge <- PercentageFeatureSet(B1_US_merge, "^HB[^(P)]", col.name = "percent_hb")

B1_US_merge<-run_rna(B1_US_merge)

p1h<-DimPlot(object = B1_US_merge, reduction = "pca", pt.size = .1, group.by = "Lane")


B1_US_merge@meta.data$Lane<-as.factor(B1_US_merge@meta.data$Lane)

B1_US_merge<-RunHarmony(B1_US_merge,group.by.vars = "Lane")

p2h<-DimPlot(object = B1_US_merge, reduction = "pca", pt.size = .1, group.by = "Lane")
library(cowplot)

pdf("harmony.pdf")
plot_grid(p1h,p2h)
dev.off()


B1_US_merge_o<-B1_US_merge

 
B1_US_merge<-run_cite(B1_US_merge)
B1_US_merge<-run_wnn(B1_US_merge)


B1_US_merge<-run_rna_harmony(B1_US_merge)
B1_US_merge<-run_cite_harmony(B1_US_merge)
B1_US_merge<-run_wnn_harmony(B1_US_merge)


seur<-B1_US_merge

library(SeuratDisk)
reference <- LoadH5Seurat("./pbmc_multimodal.h5seurat")


library(SeuratData)
InstallData('pbmc3k')

###SUPER IMPORTANT
DefaultAssay(seur)<-"SCT"

anchors <- FindTransferAnchors(
  reference = reference,
  #query = seur,
  query = seur,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)


seur <- MapQuery(
  anchorset = anchors,
  query = seur,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)




seur$mono1<-monaco_ann1(seur)

seur$mono2<-monaco_ann2(seur)

seur$atlas<-atlas_ann(seur)



#saveRDS(seur, "306_2_seur_9_16.RDS")
#saveRDS(seur, "306_2_seur_18_21.RDS")
#saveRDS(seur, "306_2_seur_22_24.RDS")

saveRDS(seur, output)
