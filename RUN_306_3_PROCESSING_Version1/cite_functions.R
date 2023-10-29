## 27 Oct 2023 Version 0.01



#' Read h5 Files
read_h5_files <- function(lanes, rawdir) {
  B1_US_data <- list()
  for(i in 1:length(lanes)){
    B1_US_data[[i]] = Read10X_h5(paste(rawdir, "23_306_3_", lanes[i], "_sample_filtered_feature_bc_matrix.h5", sep = ""))
  }
  return(B1_US_data)
}

#' Create Seurat Object
create_seurat_objects <- function(B1_US_data, lanes) {
  B1_US_SeuratObj <- list()
  for(i in 1:length(lanes)){
    B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = (B1_US_data[[i]]$'Gene Expression'), assay = "RNA", min.feature = 0)
    B1_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:142,])
    B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17), lanes[[i]], sep = ""))
    B1_US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_US_SeuratObj[[i]])))
    B1_US_SeuratObj[[i]]$Lane  <- rep(lanes[[i]], length(colnames(B1_US_SeuratObj[[i]])))
  }
  names(B1_US_data) = names(B1_US_SeuratObj) = lanes
  return(B1_US_SeuratObj)
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
read_demux_files <- function(lanes, demuxdir, B1_US_merge) {
  B1_US_demuxbestList <- list()
  for(i in 1:length(lanes)){
    B1_US_demuxbestList[[i]] = read.table(paste(demuxdir, "S2_23_306_3_",lanes[i],".best", sep = ""), sep = "\t", header = TRUE)
  }
  
  for(i in 1:length(lanes)){
    B1_US_demuxbestList[[i]]$"BEST.GUESS" =  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,B1_US_demuxbestList[[i]]$BEST.GUESS)
    B1_US_demuxbestList[[i]]$"NEXT.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,B1_US_demuxbestList[[i]]$NEXT.GUESS)
    B1_US_demuxbestList[[i]]$"SNG.BEST.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,B1_US_demuxbestList[[i]]$SNG.BEST.GUESS)
    B1_US_demuxbestList[[i]]$"SNG.NEXT.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,B1_US_demuxbestList[[i]]$SNG.NEXT.GUESS)
    B1_US_demuxbestList[[i]]$"DBL.BEST.GUESS"=  gsub("/Aligned.*|outdir/Manthiram_|_output|_S39|_S45|_S44" , '',perl=TRUE,B1_US_demuxbestList[[i]]$DBL.BEST.GUESS)
  }
  
  for(i in 1:length(lanes)){
    B1_US_demuxbestList[[i]]$NewBarcode = paste(substr(B1_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17), lanes[[i]], sep = "")
  }
  
  B1_US_demuxbestdf <- plyr::ldply(B1_US_demuxbestList, data.frame)
  rownames(B1_US_demuxbestdf) <- B1_US_demuxbestdf$NewBarcode
  B1_US_merge <- subset(B1_US_merge, cells =  B1_US_demuxbestdf$NewBarcode)
  B1_US_merge <- AddMetaData(B1_US_merge, metadata = B1_US_demuxbestdf[colnames(B1_US_merge),])
  
  return(B1_US_merge)
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
filter_cells_for_antibody_clumps <- function(B1_US_merge) {
  B1_US_merge <- subset(B1_US_merge, subset = nCount_CITE < quantile(B1_US_merge$nCount_CITE, probs = c(0.995)))
  return(B1_US_merge)
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


DSB_Normalize_and_Save <- function(B1_US_merge) {
  isotype.control.name.vec <- c(rownames(B1_US_merge[["CITE"]][c(61,63,66,117,115,97,137,83 ),]))
  # Subset data
  SNG.US <- subset(B1_US_merge, subset = DROPLET.TYPE == "SNG")
  NEG.US <- subset(B1_US_merge, subset = DROPLET.TYPE == "AMB")
  
  # Save subsets
  saveRDS(SNG.US, "SNG.US.RDS")
  saveRDS(NEG.US, "NEG.US.RDS")
  
  # Perform DSB normalization
  norm_result <- DSBNormalizeProtein(
    cell_protein_matrix = as.matrix(GetAssayData(SNG.US[["CITE"]], slot = "counts")),
    empty_drop_matrix = as.matrix(GetAssayData(NEG.US[["CITE"]], slot = "counts")),
    define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = FALSE,
    isotype.control.name.vec = isotype.control.name.vec
  )
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
           plot = DimPlot(object = B1_US_merge, reduction = "rna.umap", pt.size = .1, split.by = "BEST.GUESS")) + NoLegend()
  } else {
   
   
    ggsave(paste0(save_name,"_harmony.pdf"), 
           plot = DimPlot(object = B1_US_merge, reduction = "rna.umap.harmony", pt.size = .1, split.by = "BEST.GUESS")) + NoLegend()
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
