#' 
#' Create a pie chart that illustrates the percentage of the lineages belong to CL and tumor k nearest neighbors samples for the metasample
#'
#' @import dplyr
#' @importFrom SNFtool dist2
#' @importFrom reshape2 melt
#' @import ggplot2
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param selected_samples vector of samples that will compone the metasample
#' @param n number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @return a pie chart illustrating the distribution of the CL and TCGA k nearest neighbors lineages to the metasample 
#' @export

get_piechart_both_2 <- function(combined_mat, selected_samples, n, ann) {
  
  x_1 <- combined_mat[which(rownames(combined_mat) %in% selected_samples),]
  x_2 <- as.data.frame(colMeans(x_1))
  x_3 <- t(x_2)
  rownames(x_3) <- 'metasample'
  mat_metasample <- as.matrix(x_3)
  
  dist_metasample <- SNFtool::dist2(mat_metasample, combined_mat)
  dist_metasmaple_1 <- reshape2::melt(dist_metasample)
  colnames(dist_metasmaple_1)[c(1,2,3)] <- c('metasample','sampleID','dist')
  
  dist_top_n <- dist_metasmaple_1 %>% 
    group_by(metasample) %>%                  
    arrange(dist) %>%  
    slice_head(n = n + length(selected_samples)) %>% 
    ungroup()
  
  got_sample <- c(unique(as.character(dist_top_n$sampleID)), unique(selected_samples))
  got_sample <- got_sample[!got_sample %in% selected_samples]

  dist_df <- as.data.frame(got_sample)
  colnames(dist_df)[1] <- 'sampleID'
  
  dist_top25_1 <- left_join(dist_df, ann[,c(1,2)], by = 'sampleID')
  colnames(dist_top25_1)[c(1,2)] <- c('sample_1','lineage_tcga')
  
  dist_top25_2 <- dist_top25_1 %>% mutate('sampleID' = rep(selected_samples[1], length(dist_top25_1$sample_1))) %>% 
    left_join(., ann[,c(1,2)], by = 'sampleID')
  
  colnames(dist_top25_2)[c(3,4)] <- c('sample_2','lineage_ccle')
  
  dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, lineage_tcga, lineage_ccle)
  
  dist_top25_4 <- dist_top25_3 %>% 
    select(lineage_tcga, lineage_ccle) %>%
    table() %>% as.data.frame()
  
  dist_top25_4 <- dist_top25_4 %>%
    mutate(percentage = Freq / sum(Freq) * 100,
           label = paste0(round(percentage, 1), "%"))
  
  colnames(dist_top25_4)[1] <- 'lineage'
  
  dist_top25_4 <- dist_top25_4 %>%
    filter(Freq > 0)
  
  y <-  ggplot(dist_top25_4, aes(x = "", y = Freq, fill = lineage)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Spectral") + 
    theme_void() +      
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + # personalizza il titolo
    geom_text(aes(label = label), 
              position = position_stack(vjust = 0.5), 
              size = 5) +
    labs(title = paste("Distribution of Lineage TCGA&CCLE for metasample"))
  
  return(y)
  
}