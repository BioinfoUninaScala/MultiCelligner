#' 
#' Method to take in a Seurat object and run default Seurat clustering algorithm for mutational signature data
#' 
#' Adapted from: https://github.com/broadinstitute/celligner/blob/d9c9246f8a1b6885d07f2f28bbdca24253e57cf1/R/Celligner_methods.R
#' 
#' @import Seurat
#' @import SeuratObject
#' @param seu_obj Seurat object which contain mutational signature data
#' @return a Seurat object
#' @export

signature_cluster_data <- function (seu_obj) 
  
{
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = "pca", 
                                   dims = 1:25, k.param = 20, force.recalc = TRUE, 
                                   verbose = FALSE)
  seu_obj %<>% Seurat::FindClusters(reduction = "pca", resolution = celligner_global$mod_clust_res)
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  return(seu_obj)
}