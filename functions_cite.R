





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
  seur_obj <- RunPCA(seur_obj,reduction.name = 'pca_h')
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
seur_obj <- RunPCA(seur_obj,reduction.name = 'apca')
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
