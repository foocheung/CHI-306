## 27 Jan



generate_summary_data <- function(metadata, new_dir_name) {
  summary_data <- metadata %>%
    group_by(Lane) %>%
    summarise(mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
              mean_nCount_CITE = mean(nCount_CITE, na.rm = TRUE),
              mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
              mean_percent.mito = mean(percent.mito, na.rm = TRUE))
  
  write.csv(summary_data, file = paste(new_dir_name, "/summary_data.csv", sep = ""))
}

generate_sng_data <- function(metadata, new_dir_name) {
  sng_data <- metadata %>%
    group_by(Lane, DROPLET.TYPE) %>%
    summarise(n = n())
  
  write.csv(sng_data, file = paste(new_dir_name, "/sng_data.csv", sep = ""))
}






#' Read h5 Files
read_h5_files <- function(lanes, rawdir,prefix_dir) {
  B1_US_data <- list()
  for(i in 1:length(lanes)){
    B1_US_data[[i]] = Read10X_h5(paste(rawdir, prefix_dir, lanes[i], "_sample_filtered_feature_bc_matrix.h5", sep = ""))
  }
  return(B1_US_data)
}

#' Create Seurat Object
create_seurat_objects <- function(B1_US_data, lanes) {
  Celltype <- rep(c("Sneg","Sneg","Sneg","Spos","CD4","CD4","CD8","CD8"),3)
#  B1_US_data <- list()

   for(i in 1:length(lanes)){
     rownames(B1_US_data[[i]]$`Antibody Capture`) = gsub("\\.1$","", rownames(B1_US_data[[i]]$`Antibody Capture`)) # this regex matches all ".1" at end of string, replaces with nothing
   }
  
  SeuratObj<- list() 
for(i in 1:length(lanes)){
 
  SeuratObj[[i]] <- CreateSeuratObject(counts = B1_US_data[[i]]$'Gene Expression', assay = "RNA",min.feature = 20)

  SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:24,colnames(SeuratObj[[i]])])


  #B1_US_data[[i]][["HTO"]] <- CreateAssayObject(counts
  # = Run202308_data[[i]]$`Custom`[1:2,colnames(B1_US_data[[i]])])
  SeuratObj[[i]] <- RenameCells(SeuratObj[[i]],
                                new.names = paste(substr(colnames(SeuratObj[[i]]),
                                                         start = 1, stop = 17),lanes[[i]],"_306_4", sep = ""))
  SeuratObj[[i]]$Batch  <- rep("306_4",
                               length(colnames(SeuratObj[[i]])))
  SeuratObj[[i]]$Type  <- rep("sorted",
                              length(colnames(SeuratObj[[i]])))
  SeuratObj[[i]]$Lane <- rep(lanes[i],length(colnames(SeuratObj[[i]])))
  names(SeuratObj)[i] <- Celltype[i]
  SeuratObj[[i]]$CellType <-rep(Celltype[i],length(colnames(SeuratObj[[i]])))

 # B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = (B1_US_data[[i]]$'Gene Expression'), assay = "RNA", min.feature = 0)
 # B1_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:24,])
#  B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17), lanes[[i]], sep = ""))
#  B1_US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_US_SeuratObj[[i]])))
#  B1_US_SeuratObj[[i]]$Lane  <- rep(lanes[[i]], length(colnames(B1_US_SeuratObj[[i]])))
  
  
  
  
}
names(B1_US_data) = names(SeuratObj) = lanes
return(SeuratObj)
}




#' Merge Seurat Objects
merge_seurat_objects <- function(B1_US_SeuratObj, lanes) {
  if (length(lanes) == 1){
    B1_US_merge = B1_US_SeuratObj[[1]]
  } else {
    B1_US_merge = merge(B1_US_SeuratObj[[1]], B1_US_SeuratObj[2:length(lanes)])
  }
  return(B1_US_merge)
}

#' Add Metadata and Perform Filtering
add_metadata_and_filter <- function(B1_US_merge, save_name) {
  mito.genes = grep(pattern = "^MT-", x = rownames(B1_US_merge), value = TRUE)
  B1_US_merge <- AddMetaData(object = B1_US_merge, metadata = Matrix::colSums(B1_US_merge[mito.genes,])/Matrix::colSums(B1_US_merge), col.name = "percent.mito")
  B1_US_merge <- PercentageFeatureSet(B1_US_merge, "^RP[SL]", col.name = "percent_ribo")
  B1_US_merge <- PercentageFeatureSet(B1_US_merge, "^HB[^(P)]", col.name = "percent_hb")
  
  ggsave(paste0(save_name,"_prefiltering.pdf"), 
         plot=VlnPlot(B1_US_merge, c( "percent.mito","percent_ribo","percent_hb", "nFeature_RNA", "nCount_RNA"), 
                      group.by = "Lane", pt.size=0, ncol = 3)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  B1_US_merge <- B1_US_merge %>%
    filter(nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 20000 & percent.mito < 0.30)
  
  ggsave(paste0(save_name,"_postfiltering.pdf"), 
         plot =VlnPlot(B1_US_merge, c( "percent.mito","percent_ribo","percent_hb", "nFeature_RNA", "nCount_RNA"), 
                       group.by = "Lane", pt.size=0, ncol = 3)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
  
    return(B1_US_merge)
}

#' Read DEMUX Files
read_demux_files <- function(lanes, prefix_dir,demuxdir, B1_US_merge) {
  SeuratObj<-B1_US_merge
  
 # setwd("/Users/cheungf/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/R/Manthiram_Covid/306_4/RUN/DEMUX")
  
  B1_US_demuxbestList <- list()
  demuxbestList = list()
  for(i in 1:length(lanes)){
    demuxbestList[[i]] = read.table(paste(demuxdir,"S2_multi_config_",
                                          lanes[[i]],".best", sep = ""), sep = "\t", header = TRUE)
  }
  names(demuxbestList) <- lanes
  
  
  ##NOT GOOD TO CLEAN HERE FOR DOUBLETS REMOVES HITS
  # for(i in 1:length(Lanes)){
  #   demuxbestList[[i]]$"BEST.GUESS" =  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,demuxbestList[[i]]$BEST.GUESS)
  #   demuxbestList[[i]]$"NEXT.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,demuxbestList[[i]]$NEXT.GUESS)
  #   demuxbestList[[i]]$"SNG.BEST.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,demuxbestList[[i]]$SNG.BEST.GUESS)
  #   demuxbestList[[i]]$"SNG.NEXT.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,demuxbestList[[i]]$SNG.NEXT.GUESS)
  #   demuxbestList[[i]]$"DBL.BEST.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,demuxbestList[[i]]$DBL.BEST.GUESS)
  # }
  
  for(i in 1:length(lanes)){
    demuxbestList[[i]]$NewBarcode = paste(substr(demuxbestList[[i]]$BARCODE, 
                                                 start = 1, stop = 17),
                                          lanes[[i]],"_306_4", sep = "")
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
  
  return(SeuratObj)
}

#' Summarize Lane and Subject Data
summarize_data <- function(B1_US_merge) {
  lane_summ_dat <- B1_US_merge %>% 
    group_by(Lane, DROPLET.TYPE) %>%
    summarize(n = n())
 
  
   subj_summ_dat <- B1_US_merge %>%
    filter(DROPLET.TYPE == "SNG") %>%
    group_by(Lane, SNG.BEST.GUESS) %>%
    summarize(n = n())
  
  return(list(lane_summ_dat = lane_summ_dat, subj_summ_dat = subj_summ_dat))
}

#' Create Plots
create_plots <- function(lane_summ_dat, subj_summ_dat, save_name) {
  plot_path2 <- paste0(save_name,"_qc_plots_count_cells_filter_singlecell.pdf")
  pdf(plot_path2)
  p1 <- ggplot(lane_summ_dat, aes(x = Lane, y = n)) +
    geom_col(aes(fill= DROPLET.TYPE), position = "dodge") +
    coord_flip()
  p2 <- ggplot(subj_summ_dat, aes(y = Lane, x = SNG.BEST.GUESS)) +
    geom_tile(aes(fill= log10(n))) +
    geom_text(aes(label = round(log10(n), 2))) +
    scale_fill_viridis_c() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  print(p1)
  print(p2)
  dev.off()
}

#' Filter Cells for Antibody Clumps
filter_cells_for_antibody_clumps <- function(ObjList) {
#  B1_US_merge <- subset(B1_US_merge, subset = nCount_CITE < quantile(B1_US_merge$nCount_CITE, probs = c(0.995)))
  
  for(i in 1:length(ObjList)){
    quantile(ObjList[[i]]$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
    ObjList[[i]] <- subset(ObjList[[i]], 
                           subset = nCount_CITE < quantile(ObjList[[i]]$nCount_CITE, probs = c(0.995)))
  }
  
  
  return(ObjList)
}

#' Add Isotype Metadata
add_isotype_metadata <- function(B1_US_merge,isotype.control.name.vec ) {
 
  B1_US_merge <- AddMetaData(object = B1_US_merge, metadata = colMeans(as.matrix(GetAssayData(B1_US_merge, assay = "CITE", slot = "counts")[isotype.control.name.vec,])), col.name = "isotype.mean")
  
  B1_US_merge = subset(B1_US_merge, subset = isotype.mean < quantile(B1_US_merge$isotype.mean, probs=0.99) )
  return(B1_US_merge)
}

#' Filter Out CITEseq Outliers
filter_out_citeseq_outliers <- function(B1_US_merge) {
  B1_US_merge <- AddMetaData(object = B1_US_merge, metadata = log1p(B1_US_merge$nCount_CITE), col.name = "nCount_CITE.log1p")
  return(B1_US_merge)
}


DSB_Normalize <- function(ObjList) {
  DSB.ADT.list = list()
  
  ObjList_neg = list()
  ObjList_SNG = list()
  for(i in 1:length(ObjList)){
    ObjList_neg[[i]] <- subset(ObjList[[i]], DROPLET.TYPE == "AMB"|nFeature_RNA < 100)
    # ObjList_neg[[i]] <- subset(ObjList[[i]], cells = sample(Cells(ObjList_neg[[i]]),1000))
    ObjList_SNG[[i]] <- subset(ObjList[[i]], nFeature_RNA > 100 & DROPLET.TYPE == "SNG")
  }
  
  isotype.control.name.vec = c("MouseIgG2a", 
                               "MouseIgG2b", 
                               "MouseIgG1")
  ObjList_filt <- list()
  for(i in 1:length(ObjList)){
    DSB.ADT.list[[i]] = DSBNormalizeProtein(cell_protein_matrix = as.matrix(GetAssayData(ObjList_SNG[[i]][["CITE"]], slot = "counts")),
                                            empty_drop_matrix = as.matrix(GetAssayData(ObjList_neg[[i]][["CITE"]], slot = "counts")), 
                                            denoise.counts = TRUE, isotype.control.name.vec = isotype.control.name.vec)
    ObjList_filt[[i]] <- SetAssayData(ObjList_SNG[[i]], assay = "CITE", new.data = DSB.ADT.list[[i]], slot = "data")
  }
  
  return(ObjList_filt)
  }


filter<- function(ObjList_filt) {
for(i in 1:length(ObjList_filt)){
  ObjList_filt[[i]][["percent.mito"]] <- PercentageFeatureSet(object = ObjList_filt[[i]], 
                                                              pattern = "^MT-")
}
  
  for(i in 1:length(ObjList_filt)){
    ObjList_filt[[i]] = subset(ObjList_filt[[i]], subset = nFeature_RNA > 200)
    ObjList_filt[[i]] = subset(ObjList_filt[[i]], nFeature_RNA < 5000)
    ObjList_filt[[i]] = subset(ObjList_filt[[i]], nCount_RNA < 30000)
    ObjList_filt[[i]] = subset(ObjList_filt[[i]], percent.mito < 30)
    ObjList_filt[[i]] = subset(ObjList_filt[[i]], nCount_CITE < 20000)
  }
  
  return(ObjList_filt)
  }


#' Normalize, Find Variable Features, and Scale Data
normalize_and_scale_data <- function(B1_US_merge) {
  B1_US_merge <- NormalizeData(B1_US_merge)
  B1_US_merge <- FindVariableFeatures(B1_US_merge)
  B1_US_merge <- ScaleData(B1_US_merge)
  return(B1_US_merge)
}

#' Run CITE
run_cite_analysis <- function(B1_US_merge) {
  B1_US_merge <- run_cite(B1_US_merge) 
  return(B1_US_merge)
}

#' Run WNN
run_wnn_analysis <- function(B1_US_merge) {
  B1_US_merge <- run_wnn(B1_US_merge) 
  return(B1_US_merge)
}


#' Run RNA
run_rna_analysis <- function(B1_US_merge) {
  B1_US_merge <- run_rna(B1_US_merge) 
  return(B1_US_merge)
}


#' Run RNA harmony
#run_rna_harmony_analysis <- function(B1_US_merge) {

#  B1_US_merge <- run_rna_harmony(B1_US_merge) 
#  return(B1_US_merge)
#}

#' Generate UMAP and Harmony Plots
generate_umap_and_harmony_plots <- function(B1_US_merge, save_name, lanes) {
  B1_US_merge@meta.data$Lane <- as.factor(B1_US_merge@meta.data$Lane)
  if (length(lanes) == 1){
    ggsave(paste0(save_name,"_umap.pdf"), 
           plot = DimPlot(object = B1_US_merge,  reduction = "rna.umap.harmony", pt.size = .1, split.by = "BEST.GUESS",ncol=1, label = TRUE, repel = TRUE) & NoLegend() )
  } else {
   
   
    ggsave(paste0(save_name,"_harmony.pdf"), 
           plot = DimPlot(object = B1_US_merge,  reduction = "rna.umap.harmony", pt.size = .1, split.by = "BEST.GUESS",ncol=1, label = TRUE, repel = TRUE) & NoLegend() )
  }
}

#' Save Seurat Object as RDS
save_seurat_as_rds <- function(B1_US_merge, save_name) {
  save_name <- paste0(save_name, ".RDS")
  saveRDS(B1_US_merge, save_name)
}







run_rna <- function(seur_obj) {
 # seed(1234)
  DefaultAssay(seur_obj) <- 'RNA'
  seur_obj <- SCTransform(seur_obj)
  seur_obj <- RunPCA(seur_obj)
  seur_obj <- FindNeighbors(seur_obj, dims = 1:30)
  seur_obj <- FindClusters(seur_obj, resolution = 0.8, algorithm=3, verbose = FALSE)
  
  seur_obj <- RunUMAP(seur_obj, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                 reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

}

run_rna_harmony <- function(seur_obj) {
  # seed(1234)
  DefaultAssay(seur_obj) <- 'RNA'
  seur_obj <- SCTransform(seur_obj)
  seur_obj <- RunPCA(seur_obj,reduction.name = 'pca_h',reduction.key = 'pca_h')
  seur_obj <- FindNeighbors(seur_obj,"harmony", dims = 1:30)
  seur_obj <- FindClusters(seur_obj, resolution = 0.8, algorithm=3, verbose = FALSE)
  
  seur_obj <- RunUMAP(seur_obj, reduction = 'harmony', dims = 1:30, assay = 'RNA', 
                      reduction.name = 'rna.umap.harmony', reduction.key = 'rnaUMAP.harmony_')
  
}

run_cite_harmony <-function(seur_obj) {
  #seed(1234)
  DefaultAssay(seur_obj) <- 'CITE'
  VariableFeatures(seur_obj) <- rownames(seur_obj[["CITE"]])
  
  seur_obj <- ScaleData(seur_obj)
  seur_obj <- RunPCA(seur_obj,reduction.name = 'apca_h')
  seur_obj <- FindNeighbors(seur_obj, "harmony",dims = 1:min(length(rownames(seur_obj[['CITE']]))-1, 20), reduction = "apca")
  
  seur_obj <- FindClusters(seur_obj, graph.name = "CITE_snn", algorithm = 3, verbose = FALSE)
  seur_obj <- RunUMAP(seur_obj, reduction = 'harmony', dims = 1:30, assay = 'ADT', 
                      reduction.name = 'adt.umap.harmony', reduction.key = 'adtUMAP.harmony_')
  
}

run_cite <-function(seur_obj) {
  #seed(1234)
DefaultAssay(seur_obj) <- 'CITE'
VariableFeatures(seur_obj) <- rownames(seur_obj[["CITE"]])

seur_obj <- ScaleData(seur_obj)
seur_obj <- RunPCA(seur_obj,reduction.name = 'apca',reduction.key = 'apca_')
seur_obj <- FindNeighbors(seur_obj, dims = 1:min(length(rownames(seur_obj[['CITE']]))-1, 20), reduction = "apca")

seur_obj <- FindClusters(seur_obj, graph.name = "CITE_snn", algorithm = 3, verbose = FALSE)
seur_obj <- RunUMAP(seur_obj, reduction = 'apca', dims = 1:30, assay = 'ADT', 
               reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

}





run_wnn_harmony <-function(seur_obj) {
  
  #  seed(1234)
  seur_obj <- FindMultiModalNeighbors(
    seur_obj, reduction.list = list("pca_h", "apca_h"), 
    dims.list = list(1:20, 1:20), modality.weight.name = c("intRNA.weight", "intADT.weight"))
  
  seur_obj <- FindClusters(seur_obj, graph.name = "wsnn", verbose = FALSE)
  seur_obj <- RunUMAP(seur_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", reduction.key = "wnnUMAP.harmony_")
  
  
  
}

run_wnn <-function(seur_obj) {
 
#  seed(1234)
  seur_obj <- FindMultiModalNeighbors(
    seur_obj, reduction.list = list("pca", "apca"), 
    dims.list = list(1:20, 1:20), modality.weight.name = c("intRNA.weight", "intADT.weight"))
  
  seur_obj <- FindClusters(seur_obj, graph.name = "wsnn", verbose = FALSE)
  seur_obj <- RunUMAP(seur_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  
  
  
}




atlas_ann <-function(seur_obj) {



primary.ref <- celldex::HumanPrimaryCellAtlasData()

seur_obj <- as.SingleCellExperiment(seur_obj)
primary.main <- SingleR(test = seur_obj,assay.type.test = 1,ref = primary.ref,labels = primary.ref$label.main)

seur_obj$primary.main <- primary.main$pruned.labels
}



monaco_ann1 <-function(seur_obj) {

monaco.ref <- celldex::MonacoImmuneData()
seur_obj<- as.SingleCellExperiment(seur_obj)

monaco.main <- SingleR(test = seur_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)

seur_obj$monaco.main <- monaco.main$pruned.labels


}




monaco_ann2 <-function(seur_obj) {
  
  monaco.ref <- celldex::MonacoImmuneData()
  seur_obj<- as.SingleCellExperiment(seur_obj)
  
    monaco.fine <- SingleR(test = seur_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
  
  
  seur_obj$monaco.fine <- monaco.fine$pruned.labels
  
}



runTransferLearning <- function(reference_file, seur1) {
  reference <- LoadH5Seurat(reference_file)
  
  anchors <- FindTransferAnchors(
    reference = reference,
    query = seur1,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  
  seur1 <- MapQuery(
    anchorset = anchors,
    query = seur1,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  
  return(seur1)
}




performMitoAnalysis <- function(seur1) {
  mito.genes <- grep(pattern = "^MT-", x = rownames(seur1), value = TRUE)
  seur1 <- AddMetaData(object = seur1, metadata = Matrix::colSums(seur1[mito.genes, ]) / Matrix::colSums(seur1), col.name = "percent.mito")
  return(seur1)
}
