#' 
#' Function to calculate in parallel the Euclidean distance for each TCGA sample from each CL sample for MOFA2 matrix
#' 
#' @import foreach
#' @import doParallel
#' @import stats
#' @param MOFA_mat multiomics matrix of MOFA2 factor
#' @return a list of a list: the euclidean distance for each TCGA sample from each CL sample for MOFA2 matrix
#' @export

get_MOFA_fastdist_eu <- function(MOFA_mat) {
  
  cl <- makeCluster(50)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:nrow(MOFA_mat[grepl("TCGA", rownames(MOFA_mat)),]), .packages = 'stats') %dopar% {
    
    coppie_distanze <- list() 
    
    campione_2 <- rownames(MOFA_mat[grepl("TCGA", rownames(MOFA_mat)),])[i]
    
    for (h in 1:nrow(MOFA_mat[grepl("ACH-00", rownames(MOFA_mat)),])) {
      
      campione_1 <- rownames(MOFA_mat[grepl("ACH-00", rownames(MOFA_mat)),])[h]
      
      dist_1 <- dist(rbind(MOFA_mat[campione_1, ], MOFA_mat[campione_2, ]), method = "euclidean")
      
      coppie_distanze[[h]] <- list(CCLEsample = campione_1, TCGAsample = campione_2, dist_eu = dist_1)
    }
    return(coppie_distanze) 
  }
  
  stopCluster(cl)
  
  return(result)
}
