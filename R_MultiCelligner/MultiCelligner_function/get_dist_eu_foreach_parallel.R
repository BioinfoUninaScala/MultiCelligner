#' 
#' Function to calculate in parallel the Euclidean distance for each TCGA sample form each CL sample
#' 
#' @import foreach
#' @import doParallel
#' @import stats
#' @param mnn_res TCGA matrix corrected by MNN
#' @param CCLE_cor CCLE matrix with the frist 4 cPCs regressed out
#' @param core number of core
#' @return a list of a list: the euclidean distance for each TCGA sample from each CL sample
#' @export

get_fastdist_eu <- function(mnn_res, CCLE_cor, cores) {
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:nrow(mnn_res$corrected), .packages = 'stats') %dopar% {
    
    coppie_distanze <- list() 
    
    campione_2 <- rownames(mnn_res$corrected)[i]
    
    for (h in 1:nrow(CCLE_cor)) {
      
      campione_1 <- rownames(CCLE_cor)[h]
      
      dist_1 <- dist(rbind(CCLE_cor[campione_1, ], mnn_res$corrected[campione_2, ]), method = "euclidean")
      
      coppie_distanze[[h]] <- list(CCLEsample = campione_1, TCGAsample = campione_2, dist_eu = dist_1)
    }
    return(coppie_distanze) 
  }
  
  stopCluster(cl)
  
  return(result)
}

