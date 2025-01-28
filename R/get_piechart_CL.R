#' 
#' Create a pie chart that illustrates the percentage of the lineages belong to CL and tumor k nearest neighbors samples for the query sample
#'
#' @import dplyr
#' @import BiocNeighbors
#' @import ggplot2
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample samples of TCGA or CCLE choosed by the user 
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param BNindex a BiocNeighborIndex object containing precomputed index information
#' @param sample_order index of the original combined matrix
#' @return a pie chart illustrating the distribution of the CL k nearest neighbors lineages to the query  
#' @export

get_piechart_CL <- function(combined_mat, input_sample, k, ann, BNindex, sample_order) {
  
  query <- matrix(combined_mat[input_sample, ], nrow = 1)
  rownames(query)[1] <- rownames(combined_mat)[which(rownames(combined_mat) %in% input_sample)]
  out <- queryKNN(BNINDEX = BNindex,query = query, k = k + 100)
  
  combined_mat_1 <- combined_mat[sample_order,]
  sample_neighbors <- rownames(combined_mat_1[out$index,])
  sample_neighbors <- sample_neighbors[!sample_neighbors %in% input_sample]
  sample_neighbors_1 <- sample_neighbors[grepl('ACH-00', sample_neighbors)]
  sample_neighbors_2 <- sample_neighbors_1[1:k]
  
  ### devi fare un controllo qui se sono ordinati in base alla distanza:
  # se non lo sono, devi fare un left_join e poi arrange!!!!
  
  top_k_tumors_1 <- as.character(sample_neighbors_2)
  
  dist_df <- as.data.frame(top_k_tumors_1)
  colnames(dist_df)[1] <- 'sampleID'
  
  dist_top25_1 <- left_join(dist_df, ann[,c(1,2)], by = 'sampleID')
  colnames(dist_top25_1)[c(1,2)] <- c('sample_1','lineage_ccle')
  
  dist_top25_2 <- dist_top25_1 %>% mutate('sampleID' = rep(input_sample, length(dist_top25_1$sample_1))) %>% 
    left_join(., ann[,c(1,2)], by = 'sampleID')
  
  colnames(dist_top25_2)[c(3,4)] <- c('sample_2','lineage_tcga')
  
  dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, lineage_ccle, lineage_tcga)
  
  dist_top25_4 <- dist_top25_3 %>% 
    select(lineage_tcga, lineage_ccle) %>%
    table() %>% as.data.frame()
  
  dist_top25_4 <- dist_top25_4 %>%
    mutate(percentage = Freq / sum(Freq) * 100,
           label = paste0(round(percentage, 1), "%"))
  
  dist_top25_4 <- dist_top25_4 %>%
    dplyr::filter(Freq > 0)
  
  y <-  ggplot(dist_top25_4, aes(x = "", y = Freq, fill = lineage_ccle)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Spectral") + 
    theme_void() +      
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + # personalizza il titolo
    geom_text(aes(label = label), 
              position = position_stack(vjust = 0.5), 
              size = 5) +
    labs(title = paste("Distribution of Lineage CCLE for", 
                       unique(rownames(query)[1])))
  
  return(y)
  
}