library(foreach)
library(doParallel)

# prima avevi mnn_res$corrected, ora dato che devi dare solo la matrice combined_mat levi il subsetting con $

get_fastdist_eu_vexp <- function(mnn_res, CCLE_cor) {
  
  cl <- makeCluster(35)
  registerDoParallel(cl)
  
  result <- foreach(i = 1:nrow(mnn_res), .packages = 'stats') %dopar% {
    
    coppie_distanze <- list() # la lista dove si salvano i risultati deve essere all'interno del foreach
    
    campione_2 <- rownames(mnn_res)[i]
    
    for (h in 1:nrow(CCLE_cor)) {
      
      campione_1 <- rownames(CCLE_cor)[h]
      
      dist_1 <- dist(rbind(CCLE_cor[campione_1, ], mnn_res[campione_2, ]), method = "euclidean")
      
      coppie_distanze[[h]] <- list(CCLEsample = campione_1, TCGAsample = campione_2, dist_eu = dist_1)
    }
    return(coppie_distanze) # e deve essere restituita alla fine del for innestato più esterno !
  }
  
  stopCluster(cl)
  
  return(result)
}