
############################################################################# fun for calculate dist_eu:
##################################### this fun calculate the dist_eu between all CCLE sample and all mnn_res$corrected
####### insted of calculate just the CCLE that the MNN use for the alignment

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

############################################################################# try it in background

#CCLE_cor_meth <- readRDS("~/celligner/celligner_meth/file/get_prop_agree_file/CCLE_cor_meth.rds")
#mnn_res_meth <-  readRDS("~/celligner/celligner_meth/file/get_prop_agree_file/mnn_res_meth.rds")

#start_dist <- Sys.time()
#dist_couple <- get_dist_eu(mnn_res = mnn_res_meth, CCLE_cor = CCLE_cor_meth)
#end_dist <- Sys.time()
#time_dist_eu <- end_dist - start_dist

#saveRDS(time_dist_eu, "~/celligner/celligner_meth/script_dist&cor/get_dist_eu_v5/time_dist_eu.rds")

