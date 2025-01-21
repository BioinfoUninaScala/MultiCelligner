#' 
#' Function to calculate the Euclidean distance for each TCGA sample form each CL sample
#' 
#' @import stats
#' @param mnn_res TCGA matrix corrected by MNN
#' @param CCLE_cor CCLE matrix with the frist 4 cPCs regressed out
#' @return the euclidean distance for each TCGA sample from each CL sample
#' @export


get_dist_eu <- function(mnn_res, CCLE_cor) {
  
  coppie_distanze <- list()
  for (i in 1:nrow(mnn_res$corrected)) {
    campione_2 <- rownames(mnn_res$corrected)[i]
    
    for (h in 1:nrow(CCLE_cor)) {
      campione_1 <- rownames(CCLE_cor)[h]
      
      dist_1 <- dist(rbind(CCLE_cor[campione_1, ], mnn_res$corrected[campione_2, ]), method = "euclidean")
      
      coppie_distanze[[length(coppie_distanze) + 1]] <- list(CCLEsample = campione_1, TCGAsample = campione_2, dist_eu = dist_1)
    }
  }
  return(coppie_distanze)
}
