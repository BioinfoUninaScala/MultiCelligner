

create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  
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
    npcs = celligner_global$n_PC_dims, verbose = F
  )
  
  seu_obj %<>% Seurat::RunUMAP(
    assay = "RNA", dims = 1:celligner_global$n_PC_dims,
    reduction = "pca",
    n.neighbors = celligner_global$umap_n_neighbors,
    min.dist = celligner_global$umap_min_dist,
    metric = celligner_global$distance_metric, verbose = F
  )
  
  return(seu_obj)
}

