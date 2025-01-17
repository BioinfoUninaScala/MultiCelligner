signature_cluster_data <- function (seu_obj) 
  
{
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = "pca", 
                                   dims = 1:25, k.param = 20, force.recalc = TRUE, 
                                   verbose = FALSE)
  seu_obj %<>% Seurat::FindClusters(reduction = "pca", resolution = celligner_global$mod_clust_res)
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  return(seu_obj)
}