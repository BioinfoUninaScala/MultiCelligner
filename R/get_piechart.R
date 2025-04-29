#' 
#' Create a pie chart that illustrates the percentage of subtype belong to the k nearest neighbors 
#' 
#' @import dplyr
#' @import ggplot2
#' @import reshape
#' @import reshape2
#' @import magrittr
#' @import SNFtool
#' @import stringr
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param value string value to select to determine if the pie chart will be based on subtype percentage or lineage percentage
#' @return pie chart illustrating the distribution of k nearest neighbors subtype to the query 
#' @export

get_piechart <- function(combined_mat, input_sample = NULL, selected_samples = NULL , type, k, ann, value, dist_top_n) {
  
  top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
  dist_df <- data.frame(sampleID = top_k_tumors_1)
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','lineage_tcga')
  } else if ("Cell lines" %in% type & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Lineage_CCLE')
  } else if ("Tumors" %in% type & value %in% "lineage") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,2)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Lineage_TCGA')
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_tcga')
  } else if ("Cell lines" %in% type & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_CCLE')
  } else if ("Tumors" %in% type & value %in% "subtype") {
    dist_top25_1 <- left_join(dist_df, ann[, c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[1:2] <- c('sample_1','Subtype_TCGA')
  }
  
  if(value %in% "lineage") {
  dist_top25_2 <- dist_top25_1 %>%
    mutate(sampleID = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                          length(sample_1))) %>%
    left_join(., ann[, c(1,2)], by = 'sampleID')
  } else {
    dist_top25_2 <- dist_top25_1 %>%
      mutate(sampleID = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                            length(sample_1))) %>%
      left_join(., ann[, c(1,4)], by = 'sampleID')
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, lineage_tcga, lineage_ccle)
    colnames(dist_top25_3)[3] <- 'Lineage'
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_ccle, Lineage) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage"
  } else if ("Cell lines" %in% type & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_tcga')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Lineage_CCLE, lineage_tcga)
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_tcga, Lineage_CCLE) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage_CCLE"
  } else if ("Tumors" %in% type & value %in% "lineage") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','lineage_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Lineage_TCGA, lineage_ccle)
    dist_top25_4 <- dist_top25_3 %>%
      select(lineage_ccle, Lineage_TCGA) %>%
      table() %>% as.data.frame()
    fill_var <- "Lineage_TCGA"
  }
  
  if (all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, subtype_tcga, subtype_ccle)
    colnames(dist_top25_3)[3] <- 'Subtype'
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_ccle, Subtype) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype"
  } else if ("Cell lines" %in% type & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_tcga')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_CCLE, subtype_tcga)
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_tcga, Subtype_CCLE) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype_CCLE"
  } else if ("Tumors" %in% type & value %in% "subtype") {
    colnames(dist_top25_2)[3:4] <- c('sample_2','subtype_ccle')
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_TCGA, subtype_ccle)
    dist_top25_4 <- dist_top25_3 %>%
      select(subtype_ccle, Subtype_TCGA) %>%
      table() %>% as.data.frame()
    fill_var <- "Subtype_TCGA"
  }
  
  dist_top25_4 <- dist_top25_4 %>%
    mutate(percentage = Freq / sum(Freq) * 100,
           label = paste0(round(percentage, 1), "%")) %>%
    dplyr::filter(Freq > 0)
  
  if(all(c("Cell lines", "Tumors") %in% type) & value %in% "subtype") {
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_wrap = stringr::str_wrap(Subtype, width = 14))
    fill_var <- "Subtype_wrap"
  }
  
  if ("Cell lines" %in% type & value %in% "subtype") {
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_CCLE_wrap = stringr::str_wrap(Subtype_CCLE, width = 14))
    fill_var <- "Subtype_CCLE_wrap"
  }
  
  y <- ggplot(dist_top25_4, aes(x = "", y = Freq, fill = .data[[fill_var]])) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Spectral") +
    theme_void() +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 3) +
    labs(if(value %in% 'lineage') title = "Neighbors lineage distribution" else title = "Subtype lineage distribution") 
  
  return(y)
  
}
  