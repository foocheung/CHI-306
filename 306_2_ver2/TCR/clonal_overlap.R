library(scRepertoire)
library(Seurat)

data("contig_list")
cols_example<-as.matrix(colnames(contig_list[[1]]))
colclass_example<-sapply(contig_list[[1]], class)
# Initialize a list to store the results
immdata_10x_VDJ <- list()

# Iterate over numbers from 9 to 16
for (i in 9:16) {
  # Read the CSV file
  file_path <- paste0("306_3_TCR/2multi_config_", i, "_filtered_contig_annotations.csv")
  immdata_10x_VDJ[[i - 8]] <- read.csv(file_path, colClasses = colclass_example)
  
  # Match column names
  ind_match <- which((colnames(immdata_10x_VDJ[[i - 8]]) %in% colnames(contig_list[[i - 8]])) == TRUE)
  immdata_10x_VDJ[[i - 8]] <- immdata_10x_VDJ[[i - 8]][, c(ind_match)]
}

# Combine data
combinedVDJ <- combineTCR(immdata_10x_VDJ, samples = c("CD4_Human", "CD4_Human","CD4_Human", "CD4_Human","CD8_Human", "CD8_Human","CD8_Human", "CD8_Human"), ID = as.character(9:16))

# Modify barcodes
for (i in 1:length(combinedVDJ)) {
  combinedVDJ[[i]][["barcode"]] <- gsub("-1", "", combinedVDJ[[i]][["barcode"]])
}

##"CD4_Human_9_AAGCCGCGTCATATGC" BARCODE "CD8_Human_16_ATTATCCTCACTTATC"
################## SEE IF THERE IS AN UP TO DATE DATABASE FOR MAPPING ####################
################## DONE WITH VDJ UNTIL LATER #####################


################## PREPARE SEURAT OBJECT #####################

library(tidyverse)

seur<-readRDS("306_2_merge_with_SNP_v1.rds")


seur <-  subset(seur, subset = Lane %in% 9:16)

seur@meta.data<-seur@meta.data %>%
  mutate(
    Type = case_when(
      Lane == 9 ~  "ind_CD4",
      Lane == 10 ~ "ind_CD4",
      Lane == 11 ~ "ind_CD4",
      Lane == 12 ~ "ind_CD4",
      Lane == 13 ~ "ind_CD8",
      Lane == 14 ~ "ind_CD8",
      Lane == 15 ~ "ind_CD8",
      Lane == 16 ~ "ind_CD8",

      TRUE ~ "0"
    )
  )


seur@meta.data<-seur@meta.data %>%
  mutate(
    orig.ident = case_when(
     
      Lane == 9 ~  "Human_CD4_sc",
      Lane == 10 ~ "Human_CD4_sc",
      Lane == 11 ~ "Human_CD4_sc",
      Lane == 12 ~ "Human_CD4_sc",
      Lane == 13 ~ "Human_CD8_sc",
      Lane == 14 ~ "Human_CD8_sc",
      Lane == 15 ~ "Human_CD8_sc",
      Lane == 16 ~ "Human_CD8_sc",
      TRUE ~ "0"
    )
  )




cell.names <- rownames(seur@meta.data)

ind_group1<-grep('CD4', seur@meta.data[["orig.ident"]])
ind_group2<-grep('CD8', seur@meta.data[["orig.ident"]])

cell.names[ind_group1] <- paste0("CD4_Human_", cell.names[ind_group1])
cell.names[ind_group2] <- paste0("CD8_Human_", cell.names[ind_group2])

cell.names.matrix<-cell.names

#cell.names.matrix[grep("7",cell.names.matrix)]<-gsub("CD4_Human_","CD4Tcells_7_",cell.names.matrix[grep("7",cell.names.matrix)])
#cell.names.matrix[grep("8",cell.names.matrix)]<-gsub("CD8_Human_","CD8Tcells_8_",cell.names.matrix[grep("8",cell.names.matrix)])

cell.names.matrix[grep("9",cell.names.matrix)]<-gsub("CD4_Human_","CD4_Human_9_",cell.names.matrix[grep("9",cell.names.matrix)])
cell.names.matrix[grep("10",cell.names.matrix)]<-gsub("CD4_Human_","CD4_Human_10_",cell.names.matrix[grep("10",cell.names.matrix)])
cell.names.matrix[grep("11",cell.names.matrix)]<-gsub("CD4_Human_","CD4_Human_11_",cell.names.matrix[grep("11",cell.names.matrix)])
cell.names.matrix[grep("12",cell.names.matrix)]<-gsub("CD4_Human_","CD4_Human_12_",cell.names.matrix[grep("12",cell.names.matrix)])
cell.names.matrix[grep("13",cell.names.matrix)]<-gsub("CD8_Human_","CD8_Human_13_",cell.names.matrix[grep("13",cell.names.matrix)])
cell.names.matrix[grep("14",cell.names.matrix)]<-gsub("CD8_Human_","CD8_Human_14_",cell.names.matrix[grep("14",cell.names.matrix)])
cell.names.matrix[grep("15",cell.names.matrix)]<-gsub("CD8_Human_","CD8_Human_15_",cell.names.matrix[grep("15",cell.names.matrix)])
cell.names.matrix[grep("16",cell.names.matrix)]<-gsub("CD8_Human_","CD8_Human_16_",cell.names.matrix[grep("16",cell.names.matrix)])


cell.names.matrix[grep("-9",cell.names.matrix)]<-gsub("-9_306_2","",cell.names.matrix[grep("-9",cell.names.matrix)])
cell.names.matrix[grep("-10",cell.names.matrix)]<-gsub("-10_306_2","",cell.names.matrix[grep("-10",cell.names.matrix)])
cell.names.matrix[grep("-11",cell.names.matrix)]<-gsub("-11_306_2","",cell.names.matrix[grep("-11",cell.names.matrix)])
cell.names.matrix[grep("-12",cell.names.matrix)]<-gsub("-12_306_2","",cell.names.matrix[grep("-12",cell.names.matrix)])
cell.names.matrix[grep("-13",cell.names.matrix)]<-gsub("-13_306_2","",cell.names.matrix[grep("-13",cell.names.matrix)])
cell.names.matrix[grep("-14",cell.names.matrix)]<-gsub("-14_306_2","",cell.names.matrix[grep("-14",cell.names.matrix)])
cell.names.matrix[grep("-15",cell.names.matrix)]<-gsub("-15_306_2","",cell.names.matrix[grep("-15",cell.names.matrix)])
cell.names.matrix[grep("-16",cell.names.matrix)]<-gsub("-16_306_2","",cell.names.matrix[grep("-16",cell.names.matrix)])

seur<-RenameCells(seur, new.names = cell.names.matrix)

combo<-combineExpression(combinedVDJ, seur,cloneCall="aa", group.by = "none")

combo@meta.data$BEST.GUESS=  gsub("outdir/Manthiram_|_S.*" , '',perl=TRUE,combo@meta.data$BEST.GUESS)

combo@meta.data$Type_Sample<-paste0(combo@meta.data$Type,"_",combo@meta.data$BEST.GUESS,"_",combo$orig.ident)
combo@meta.data$Type_Group<-combo@meta.data$Type_Sample


combo@meta.data$Type_Group<-gsub(".*T101.*Human","INF",combo@meta.data$Type_Group)
combo@meta.data$Type_Group<-gsub(".*T103.*Human","INF",combo@meta.data$Type_Group)
combo@meta.data$Type_Group<-gsub(".*T126.*Human","VACC",combo@meta.data$Type_Group)
combo@meta.data$Type_Group<-gsub(".*T133.*Human","VACC",combo@meta.data$Type_Group)
combo@meta.data$Type_Group<-gsub(".*T10.*Human","CONT",combo@meta.data$Type_Group)

Human_CD4_integrated <- subset(combo, orig.ident == "Human_CD4_sc")
Human_CD8_integrated <- subset(combo, orig.ident == "Human_CD8_sc")


Human_CD8_integrated_sc <- as.SingleCellExperiment(Human_CD8_integrated)
##
#devtools::install_github("ncborcherding/scRepertoire@v1")
#library(scRepertoire)
#newList_Type_Group_CD8 <-expression2List(Human_CD8_integrated_sc, split.by =  "Type_Group")

p1<-clonalOverlap(Human_CD8_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita")


Human_CD4_integrated_sc <- as.SingleCellExperiment(Human_CD4_integrated)

p2<-clonalOverlap(Human_CD4_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita")

#clonalOverlap(Human_CD4_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita",exportTable = TRUE)
#CONT_CD4_sc INF_CD4_sc VACC_CD4_sc
#CONT_CD4_sc          NA 0.00071179 0.000604196
#INF_CD4_sc           NA         NA 0.000811190
#VACC_CD4_sc          NA         NA          NA

#clonalOverlap(Human_CD8_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita", exportTable = TRUE)
#CONT_CD8_sc  INF_CD8_sc VACC_CD8_sc
#CONT_CD8_sc          NA 5.37379e-05 0.000000000
#INF_CD8_sc           NA          NA 0.000216883
#VACC_CD8_sc          NA          NA          NA

pdf("clonalOverlap.pdf")
gridExtra::grid.arrange(p1,p2)
dev.off()

