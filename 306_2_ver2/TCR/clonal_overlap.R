setwd("/Users/cheungf/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/R/TCR_306_2")
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
  ind_match <- which((colnames(immdata_10x_VDJ[[i - 8]]) %in% colnames(contig_list[[1]])) == TRUE)
  immdata_10x_VDJ[[i - 8]] <- immdata_10x_VDJ[[i - 8]][, c(ind_match)]
}

# Combine data
combinedVDJ <- combineTCR(immdata_10x_VDJ, 
                          samples = c("CD4_Human", "CD4_Human","CD4_Human", "CD4_Human","CD8_Human", "CD8_Human","CD8_Human", "CD8_Human"), 
                          ID = as.character(9:16)) #, removeMulti = TRUE,removeNA = TRUE)
                        #,
                        #  removeNA = FALSE, 
                        #  removeMulti = FALSE, 
                        #  removeNA = TRUE, 
                        #  removeMulti = TRUE, 
                        #  filterMulti = TRUE)
                          




# Modify barcodes
for (i in 1:length(combinedVDJ)) {
  combinedVDJ[[i]][["barcode"]] <- gsub("-1", "", combinedVDJ[[i]][["barcode"]])
}

##"CD4_Human_9_AAGCCGCGTCATATGC" BARCODE "CD8_Human_16_ATTATCCTCACTTATC"
################## SEE IF THERE IS AN UP TO DATE DATABASE FOR MAPPING ####################
################## DONE WITH VDJ UNTIL LATER #####################


################## PREPARE SEURAT OBJECT #####################

library(tidyverse)

seur1<-readRDS("306_2_merge_with_SNP_v1.rds")


seur <-  subset(seur1, subset = Lane %in% 9:16)

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



######CHECK THE JOIN IS RIGHT HERE MIGHT BE WRONG!!!!!
#################
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

Human_CD4_integrated_sc <- as.SingleCellExperiment(Human_CD4_integrated)
Human_CD8_integrated_sc <- as.SingleCellExperiment(Human_CD8_integrated)
##
#devtools::install_github("ncborcherding/scRepertoire@v1")
#library(scRepertoire)
#newList_Type_Group_CD8 <-expression2List(Human_CD8_integrated_sc, split.by =  "Type_Group")

p1<-clonalOverlap(Human_CD4_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita",chain="TRB", )

p1B<-clonalOverlap(Human_CD4_integrated_sc,group.by = "BEST.GUESS", cloneCall = "aa", method = "raw",chain="TRB")

p11<-clonalOverlap(Human_CD4_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "morisita",chain="TRA", )

p1BB<-clonalOverlap(Human_CD4_integrated_sc,group.by = "BEST.GUESS", cloneCall = "aa", method = "raw",chain="TRA")


p1C<-scRepertoire::clonalQuant(Human_CD4_integrated_sc, cloneCall="gene+nt", scale = T, group.by = "BEST.GUESS")

p1D<-scRepertoire::clonalQuant(Human_CD4_integrated_sc, cloneCall="gene+nt", scale = F, group.by = "BEST.GUESS")

p1E<-clonalAbundance(Human_CD4_integrated_sc, cloneCall = "aa", scale = F, group.by = "BEST.GUESS")


p1F<-clonalCompare(Human_CD4_integrated_sc, 
              top.clones = 50, 
              samples = NULL, 
              cloneCall="aa", 
              graph = "alluvial", group.by = "BEST.GUESS", exportTable=F, chain="TRB") + theme(legend.position="none")

p1G<-clonalCompare(Human_CD4_integrated_sc, 
              top.clones = 50, 
              samples = NULL, 
              cloneCall="aa", 
              graph = "alluvial", group.by = "Type_Group", exportTable=F, chain="TRB") + theme(legend.position="none")



p1FF<-clonalCompare(Human_CD4_integrated_sc, 
                   top.clones = 50, 
                   samples = NULL, 
                   cloneCall="aa", 
                   graph = "alluvial", group.by = "BEST.GUESS", exportTable=F, chain="TRA") + theme(legend.position="none")

p1GG<-clonalCompare(Human_CD4_integrated_sc, 
                   top.clones = 50, 
                   samples = NULL, 
                   cloneCall="aa", 
                   graph = "alluvial", group.by = "Type_Group", exportTable=F, chain="TRA") + theme(legend.position="none")


##p1F<-clonalLength(Human_CD4_integrated_sc, cloneCall="aa", chain = "both", group.by = "Type_Group") 


#p2<-clonalOverlap(Human_CD8_integrated_sc,group.by = "Type_Group", cloneCall = "aa", method = "raw")#

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

#pdf("clonalOverlap.pdf")
#gridExtra::grid.arrange(p1,p2)
#dev.off()



p3<-clonalHomeostasis(Human_CD4_integrated_sc, 
                  cloneCall = "gene", group.by = "Type_Group")


p33<-clonalHomeostasis(Human_CD4_integrated_sc, 
                      cloneCall = "gene", group.by = "BEST.GUESS")


p4<-clonalProportion(Human_CD4_integrated_sc, 
                 cloneCall = "gene",  group.by = "Type_Group")

p44<-clonalProportion(Human_CD4_integrated_sc, 
                     cloneCall = "gene",  group.by = "BEST.GUESS")

p5<-percentAA(Human_CD4_integrated_sc, 
          chain = "TRB", 
          aa.length = 20, group.by = "Type_Group")

p6<-positionalEntropy(Human_CD4_integrated_sc, 
                chain = "TRB", 
                aa.length = 20, group.by = "Type_Group")


p7<-positionalProperty(Human_CD4_integrated_sc[c(1,2)], 
                chain = "TRB", 
                aa.length = 20, 
                method = "Atchley",  group.by = "Type_Group") 


p8<-vizGenes(Human_CD4_integrated_sc, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE,  group.by = "Type_Group")



#p9<-vizGenes(Human_CD4_integrated_sc, 
#         x.axis = "TRBV",
#         y.axis = "TRBJ",
#         plot = "heatmap",  
#         scale = TRUE,  group.by = "Type_Group")



df.genes <-percentGenes(Human_CD4_integrated_sc, 
                        chain = "TRB", 
                        gene = "Vgene", group.by = "Type_Sample", 
                         exportTable = TRUE)



#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
p2<-ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic() 



p9<-clonalDiversity(Human_CD4_integrated_sc, 
                cloneCall = "gene", group.by = "Type_Group")


p10<-clonalSizeDistribution(Human_CD4_integrated_sc, 
                       cloneCall = "aa", 
                       method= "ward.D2", group.by = "Type_Group")


#grid.arrange(bp,                                    # bar plot spaning two columns
#             bxp, sp,                               # box plot and scatter plot
#             ncol = 2, nrow = 2, 
#             layout_matrix = rbind(c(1,1), c(2,3)))


library(gridExtra)
library(patchwork)
pdf("allcd4.PDF", width=24,height=12)
library(cowplot)
grid.arrange(p1,p1B,p1F,p1G,
             p1C,p1D,p1E,p2,
             p3,p4,p5,
             p6,p7,p8,
             p9,p10)

dev.off()


pdf("allcd4.PDF", width=32,height=16)

par( mar=c(2,2,2,2) )

plot_grid(p1,p1B,p11,p1BB, 
          p1F,p1G,p1FF,p1GG,
          p1C,p1D,p1E,p2,
          p3,p33,p4,p44,p5,
          p6,p7,p8,
          p9,p10, 
          labels = c("clonalOverlap-TRB", "clonalOverlap-TRB", "clonalOverlap-TRA","clonalOverlap-TRA", 
                      "clonalCompare-TRB", "clonalCompare-TRB", "clonalCompare-TRA", "clonalCompare-TRA",
                     "clonalQuant","clonalQuant","clonalAbundance","PCA", "clonalHomeostasis","clonalHomeostasis", 
                     "clonalProportion", "clonalProportion", "percentAA", "positionalEntropy", "positionalProperty", 
                     "vizGenes", "clonalDiversity", "clonalSizeDistribution"),
          hjust = 0,vjust = 1, label_x = 0.05, scale=0.90)

dev.off()
          

saveRDS(combo, "combo.rds")
