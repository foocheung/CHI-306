library(docopt)
library(Seurat)
library(dplyr)
library(matrixStats)
library(tidyverse)
library(tidyseurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(SingleR)
library(SeuratData)
library(SeuratDisk)
InstallData('pbmc3k')


source("./functions_cite.R")

# Define the command-line interface
doc <- "
Usage:
  run_script.R [--input-dir=<input_dir>] [--output-dir=<output_dir>] [--threshold=<threshold>] [--reference-file=<reference_file>] [--output-rds=<output_rds>]

Options:
  --input-dir=<input_dir>        Input directory [default: ./DATA3/sample_filtered_feature_bc_matrix_10000_low/]
  --output-dir=<output_dir>      Output directory [default: ./RDS]
  --threshold=<threshold>        Threshold for HTO counts [default: 100]
  --reference-file=<reference_file>   Path to the reference file [default: ./pbmc_multimodal.h5seurat]
  --output-rds=<output_rds>      Output RDS file name [default: seur1_1000_9000_filter_HTO.rds]
"

# Parse the command-line arguments
args <- docopt(doc)


# Read data
rdata <- Read10X(args$'--input-dir')

seur1 <- CreateSeuratObject(counts = rdata$`Gene Expression`, assay = "RNA", min.feature = 0)
seur1[['HTO']] <- CreateAssayObject(counts = rdata$Custom)
seur1$lane <- rep("Lane1", ncol(seur1))

DefaultAssay(seur1) <- "RNA"

performMitoAnalysis(seur1)
# Run preprocessing steps
seur1 <- run_rna(seur1)

# Load reference dataset
reference_file <- args$'--reference-file'

runTransferLearning(reference_file, seur1)

# Subset based on HTO counts
seur1 <- subset(seur1, subset = nCount_HTO > as.numeric(args$'--threshold'))

# Save Seurat object as RDS file
output_dir <- args$'--output-dir'
output_file <- file.path(output_dir, args$'--output-rds')
saveRDS(seur1, file = output_file)
