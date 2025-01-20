#' 
#' Method to create seurat objects given Mutational signature matrix and annotation table
#' 
#' Adapted from: https://github.com/broadinstitute/celligner/blob/d9c9246f8a1b6885d07f2f28bbdca24253e57cf1/R/Celligner_methods.R
#' 
#' @import Seurat
#' @import SeuratObject
#' @import magrittr
#' @import tidyverse 
#' @param mat matrix of mutational signature data
#' @param ann annotation file: Expects column 'sampleID' which matches the rownames of mat.
#' @param type optional parameter: string specifying the data type of the current data (ex. 'tumor'), which is added to the annotation matrix.
#' @return a Seurat object with mutational signature data
#' @export

createSeuObj_signature <- function(mat, ann, type = NULL) {
  
  options(Seurat.object.assay.version = "v3")
  
  seu_obj <- Seurat::CreateSeuratObject(t(mat),
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


