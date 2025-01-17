#' 
#' Create a datatable of CL k nearest neighbors from the user query sorted by distance, with samples information: lineage and subtype
#'
#' @import dplyr
#' @import BiocNeighbors
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample samples of TCGA or CCLE choosed by the user 
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param BNindex a BiocNeighborIndex object containing precomputed index information
#' @param sample_order index of the original combined matrix
#' @return a datatable of CL samples with lineage and subtype information 
#' @export


get_neighbors_table_CL <- function(combined_mat, input_sample, k, ann, BNindex, sample_order) {
  
  query <- matrix(combined_mat[input_sample, ], nrow = 1)
  rownames(query)[1] <- rownames(combined_mat)[which(rownames(combined_mat) %in% input_sample)]
  out <- queryKNN(BNINDEX = BNindex,query = query, k = k + 100)
  
  combined_mat_1 <- combined_mat[sample_order,]
  sample_neighbors <- rownames(combined_mat_1[out$index,])
  sample_neighbors <- sample_neighbors[!sample_neighbors %in% input_sample]
  sample_neighbors_1 <- sample_neighbors[grepl('ACH-00', sample_neighbors)]
  sample_neighbors_2 <- sample_neighbors_1[1:k]
  
  top_k_tumors_1 <- as.character(sample_neighbors_2)
  
  dist_df <- as.data.frame(top_k_tumors_1)
  colnames(dist_df)[1] <- 'sampleID'
  
  top_k_tumors <- left_join(dist_df, ann[,c(1,2,3)], by = 'sampleID')
  top_k_tumors <- top_k_tumors %>% mutate('n' = c(1:nrow(top_k_tumors)))
  colnames(top_k_tumors)[c(1,2)] <- c('top_neighbors','lineage')
  top_k_tumors <- top_k_tumors %>% select(n,top_neighbors,lineage,subtype)
  
  return(top_k_tumors)
  
}