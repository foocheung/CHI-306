library("Seurat") #load Seurat 4
library("dplyr")
library("matrixStats")
library('tidyverse')
library(tidyseurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library("ggpubr")

# source("./functions_cite.R")

# Define docopt arguments
doc <- "
Usage:
  script.R [--rawdir=<rawdir>] [--B1LanesUnsorted=<lanes>]

Options:
  -h --help           Show this help message and exit.
  --rawdir=<rawdir>   Path to raw data directory.
  --B1LanesUnsorted=<lanes>  List of B1 lanes (comma separated).
"

# Load docopt and parse arguments
if (interactive()) {
  library(docopt)
  args <- docopt(doc)
  rawdir <- args$--rawdir
  B1LanesUnsorted <- as.numeric(unlist(strsplit(args$--B1LanesUnsorted, ",")))
} else {
  rawdir <- commandArgs(trailingOnly = TRUE)[1]
  B1LanesUnsorted <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
}

B1_US_data <- list()

for(i in seq_along(B1LanesUnsorted)){
  B1_US_data[[i]] <- Read10X_h5(paste(rawdir, "23_306_3_", B1LanesUnsorted[i], "_sample_filtered_feature_bc_matrix.h5", sep = ""))
}

B1_US_SeuratObj <- list()

for(i in seq_along(B1LanesUnsorted)){
  B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = (B1_US_data[[i]]$'Gene Expression'), assay = "RNA", min.feature = 0)
  B1_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:142,])
  B1_US_SeuratObj[[i]][["Custom"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Custom'[1:3,])
  B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17), B1LanesUnsorted[i], sep = ""))
  B1_US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_US_SeuratObj[[i]])))
  B1_US_SeuratObj[[i]]$Lane  <- rep(B1LanesUnsorted[i], length(colnames(B1_US_SeuratObj[[i]])))
}

names(B1_US_data) <- names(B1_US_SeuratObj) <- B1LanesUnsorted

# Assuming B1_US_merge is defined earlier
# Add the following lines if not defined:
# B1_US_merge <- Seurat::MergeSeurat(object.list = B1_US_SeuratObj, project = "B1UnSort")
# B1_US_merge <- Seurat::NormalizeData(B1_US_merge)

B1_US_merge$log10_CR2_plus_1 <- log10((B1_US_merge@assays$CITE@counts['CR2.1', ])+1)
B1_US_merge$log10_CXCR3_plus_1 <- log10((B1_US_merge@assays$CITE@counts['CXCR3.1', ])+1)
B1_US_merge$CR2 <- B1_US_merge@assays$CITE@counts['CR2.1', ]
B1_US_merge$CXCR3 <- B1_US_merge@assays$CITE@counts['CXCR3.1', ]

cxcr3_counts_gr3 <- B1_US_merge %>% filter(CXCR3 > 3) %>% select(CXCR3, Lane) %>% group_by(Lane) %>% tally() %>% select(2)
cxcr3_counts_gr10 <- B1_US_merge %>% filter(CXCR3 > 10) %>% select(CXCR3, Lane) %>% group_by(Lane) %>% tally() %>% select(2)

cr2_counts_gr3 <- B1_US_merge %>% filter(CR2 > 3) %>% select(CR2, Lane) %>% group_by(Lane) %>% tally() %>% select(2)
cr2_counts_gr10 <- B1_US_merge %>% filter(CR2 > 10) %>% select(CR2, Lane) %>% group_by(Lane) %>% tally() %>% select(2)

p1 <- ggplot(B1_US_merge, aes(x = log10_CXCR3_plus_1)) + geom_histogram()

p2 <- ggplot(B1_US_merge, aes(x = log10_CXCR3_plus_1)) + geom_histogram() + facet_wrap( .~ Lane ,switch = "x")

p3 <- ggplot(B1_US_merge, aes(x = log10_CR2_plus_1))  + geom_histogram(breaks = seq(0, max(3), 0.20))

p4 <- ggplot(B1_US_merge, aes(x = log10_CR2_plus_1))  + geom_histogram(breaks = seq(0, max(3), 0.20)) + facet_wrap( .~ Lane ,switch = "x")

ggexport(p1, p2, p3, p4, ncol = 2, 
         labels = c("All Samples Combined (9:16,17:24)", "Individuals Samples (9:16,17:24)"),
         filename = "counts.pdf")
