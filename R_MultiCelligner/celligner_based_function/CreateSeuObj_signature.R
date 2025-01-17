createSeuObj_signature <- function(exp_mat, ann, type = NULL) {
  
  options(Seurat.object.assay.version = "v3")
  
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$sampleID)
  )
  if (!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(
    assay = "RNA",
    features = rownames(Seurat::GetAssayData(seu_obj)),
    npcs = 25, verbose = F
  )
  
  seu_obj %<>% Seurat::RunUMAP(
    assay = "RNA", dims = 1:25,
    umap.method = "uwot",
    n.neighbors = 15,
    min.dist = 0.01,
    metric = "cosine", verbose = F
  )
  
  return(seu_obj)
}


