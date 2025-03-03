#' 
#' Create a pie chart that illustrates the percentage of subtype belong to the k nearest neighbors 
#' 
#' @import dplyr
#' @import SNFtool
#' @import ggplot2
#' @import reshape
#' @import reshape2
#' @import tidyverse
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @return a pie chart illustrating the distribution of k nearest neighbors subtype to the query 
#' @export

get_piechart_subtype <- function(combined_mat, input_sample = NULL, selected_samples = NULL , type, k, ann) {
  
  if(is.null(selected_samples)) {
    
    x_1 <- combined_mat[which(rownames(combined_mat) %in% input_sample),]
    x_2 <- as.matrix(x_1)
    x_3 <- t(x_2)
    rownames(x_3) <- input_sample
    
    dist <- SNFtool::dist2(x_3, combined_mat)
    dist_1 <- reshape2::melt(dist)
    colnames(dist_1)[c(1,2,3)] <- c('ref_ID', 'sampleID', 'dist')
    
    if ("Tumors" %in% type) { 
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>%
        mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
        filter(grepl('TCGA', x = sampleID) & priority == 2) 
    }
    
    if ("Cell lines" %in% type) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>%
        mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
        filter(grepl('TCGA', x = sampleID) & priority == 2) 
    }
    
    if (all(c("Cell lines", "Tumors") %in% type)) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>%
        mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
        filter(priority == 2)
    }
    
    dist_top_n <- dist_2 %>% 
      group_by(ref_ID) %>%                  
      arrange(dist) %>%  
      slice_head(n = k) %>% 
      ungroup()
    
    top_k_tumors_1 <- as.character(dist_top_n$sampleID)
    
  } else if (is.null(input_sample)) {
    
    x_1 <- combined_mat[which(rownames(combined_mat) %in% selected_samples),]
    x_2 <- as.data.frame(colMeans(x_1))
    x_3 <- t(x_2)
    rownames(x_3) <- 'metasample'
    mat_metasample <- as.matrix(x_3)
    
    dist_metasample <- SNFtool::dist2(mat_metasample, combined_mat)
    dist_1 <- reshape2::melt(dist_metasample)
    colnames(dist_1)[c(1,2,3)] <- c('ref_ID','sampleID','dist')
    
    if ("Tumors" %in% type) { 
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(grepl('TCGA', x = sampleID) & priority == 2)
    } 
    
    if ("Cell lines" %in% type) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(grepl('ACH-00', x = sampleID) & priority == 2)
    }
    
    if (all(c("Cell lines", "Tumors") %in% type)) {
      
      dist_2 <- dist_1 %>% 
        arrange(dist) %>% 
        mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
        filter(priority == 2)
    }
  }

    dist_top_n <- dist_2 %>% 
    group_by(ref_ID) %>%                  
    arrange(dist) %>%  
    slice_head(n = k) %>% 
    ungroup()
  
  top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
  
  if ("Cell lines" %in% type) {
    
    dist_df <- as.data.frame(top_k_tumors_1)
    colnames(dist_df)[1] <- 'sampleID'
    
    dist_top25_1 <- left_join(dist_df, ann[,c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[c(1,2)] <- c('sample_1','Subtype_CCLE')
    
    dist_top25_2 <- dist_top25_1 %>% mutate('sampleID' = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                                                             length(dist_top25_1$sample_1))) %>% 
      left_join(., ann[,c(1,2)], by = 'sampleID') 
    
    colnames(dist_top25_2)[c(3,4)] <- c('sample_2','Subtype_TCGA')
    
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_CCLE, Subtype_TCGA)
    
    dist_top25_4 <- dist_top25_3 %>% 
      select(Subtype_TCGA, Subtype_CCLE) %>%
      table() %>% as.data.frame()
    
    dist_top25_4 <- dist_top25_4 %>%
      mutate(percentage = Freq / sum(Freq) * 100,
             label = paste0(round(percentage, 1), "%"))
    
    dist_top25_4 <- dist_top25_4 %>%
      dplyr::filter(Freq > 0)
    
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_CCLE_wrapped = stringr::str_wrap(Subtype_CCLE, width = 14))  
    
    y <-  ggplot(dist_top25_4, aes(x = "", y = Freq, fill = Subtype_CCLE_wrapped)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      scale_fill_brewer(palette = "Spectral") + 
      theme_void() +      
      theme(axis.text.x = element_blank(),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + 
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      labs(title = "Neighbors subtype distribution")
    
  } 
  
  if ("Tumors" %in% type) {
    
    dist_df <- as.data.frame(top_k_tumors_1)
    colnames(dist_df)[1] <- 'sampleID'
    
    dist_top25_1 <- left_join(dist_df, ann[,c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[c(1,2)] <- c('sample_1','Subtype_TCGA')
    
    dist_top25_2 <- dist_top25_1 %>% mutate('sampleID' = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                                                             length(dist_top25_1$sample_1))) %>% 
      left_join(., ann[,c(1,4)], by = 'sampleID') 
    
    colnames(dist_top25_2)[c(3,4)] <- c('sample_2','Subtype_CCLE')
    
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype_TCGA, Subtype_CCLE)
    
    dist_top25_4 <- dist_top25_3 %>% 
      select(Subtype_CCLE, Subtype_TCGA) %>%
      table() %>% as.data.frame()
    
    dist_top25_4 <- dist_top25_4 %>%
      mutate(percentage = Freq / sum(Freq) * 100,
             label = paste0(round(percentage, 1), "%"))
    
    dist_top25_4 <- dist_top25_4 %>%
      dplyr::filter(Freq > 0)
    
    y <-  ggplot(dist_top25_4, aes(x = "", y = Freq, fill = Subtype_TCGA)) +
      geom_bar(width = 1, stat = "identity") +
      scale_y_continuous(expand = c(0, 0)) +
      coord_polar("y", start = 0) +
      scale_fill_brewer(palette = "Spectral") + 
      theme_void() +      
      theme(axis.text.x = element_blank(),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + 
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      labs(title = "Neighbors subtype distribution")
    
  }
  
  if (all(c("Cell lines", "Tumors") %in% type)) {
    
    dist_df <- as.data.frame(top_k_tumors_1)
    colnames(dist_df)[1] <- 'sampleID'
    
    dist_top25_1 <- left_join(dist_df, ann[,c(1,4)], by = 'sampleID')
    colnames(dist_top25_1)[c(1,2)] <- c('sample_1','Subtype')
    
    dist_top25_2 <- dist_top25_1 %>% mutate('sampleID' = rep(if (is.null(selected_samples)) input_sample else selected_samples[1], 
                                                             length(dist_top25_1$sample_1))) %>% 
                                                              left_join(., ann[,c(1,4)], by = 'sampleID') 
    
    colnames(dist_top25_2)[c(3,4)] <- c('sample_2','Subtype_CCLE')
    
    dist_top25_3 <- dist_top25_2 %>% select(sample_1, sample_2, Subtype, Subtype_CCLE)
    
    dist_top25_4 <- dist_top25_3 %>% 
      select(Subtype_CCLE, Subtype) %>%
      table() %>% as.data.frame()
    
    dist_top25_4 <- dist_top25_4 %>%
      mutate(percentage = Freq / sum(Freq) * 100,
             label = paste0(round(percentage, 1), "%"))
    
    dist_top25_4 <- dist_top25_4 %>%
      dplyr::filter(Freq > 0)
    
    dist_top25_4 <- dist_top25_4 %>%
      mutate(Subtype_CCLE_wrapped = stringr::str_wrap(Subtype_CCLE, width = 14))  
    
    y <-  ggplot(dist_top25_4, aes(x = "", y = Freq, fill = Subtype)) +
      geom_bar(width = 1, stat = "identity") +
      scale_y_continuous(expand = c(0, 0)) +
      coord_polar("y", start = 0) +
      scale_fill_brewer(palette = "Spectral") + 
      theme_void() +      
      theme(axis.text.x = element_blank(),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) + 
      geom_text(aes(label = label), 
                position = position_stack(vjust = 0.5), 
                size = 3) +
      labs(title = "Neighbors subtype distribution")
    
  }
  
  return(y)
  
}
