#' 
#' Create an interactive plot that highlighting k nearest neighbors for the choosed query
#'
#' @import plotly
#' @import crosstalk
#' @import reactable
#' @import htmltools
#' @import fontawesome
#' @import dplyr
#' @import magrittr
#' @import reshape
#' @import SNFtool
#' @import reshape2
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param input_sample single TCGA or CCLE sample chosen by the user 
#' @param selected_samples multiple TCGA or CCLE samples chosen by the user
#' @param type type of the k neighbors samples (Tumor or Cell Line)
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param query_lineage limit the neighbors searching to that lineage/s
#' @return an interactive plot that highlighting the tumor k nearest neighbors 
#' @export
#' 

find_neighbors <- function(combined_mat, reduced_mat, input_sample = NULL, selected_samples = NULL, type, k, ann, query_lineage) {

    if ("All" %in% query_lineage && length(query_lineage) > 1) {
    showNotification("Select only All or choose multiple lineage without All")
    warning("Select only All or choose multiple lineage without All")
    return(NULL)}
  
  ### when you switch omics with the sample, if the there isn't the sample in that omics will appear a notification!
  if(!is.null(input_sample)) {
    if(all(!input_sample %in% colnames(reduced_mat))) {
      showNotification("There are no tumor samples for this lineage in this omics")
      warning("There are no tumor samples for this lineage in this omics")
      return(NULL)
    }
  }
  
  ### when you switch omics with the sample, if the there isn't the sample in that omics will appear a notification!
  if(!is.null(selected_samples)) {
    if(all(!selected_samples %in% colnames(reduced_mat))) {
      showNotification("There are no tumor samples for this lineage in this omics")
      warning("There are no tumor samples for this lineage in this omics")
      return(NULL)
    }
  }
  
  if(all(query_lineage == 'All')) {
    
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
          filter(grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', x = sampleID) & priority == 2) 
      }
      
      if ("Cell lines" %in% type) {
        
        dist_2 <- dist_1 %>% 
          arrange(dist) %>%
          mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
          filter(grepl('ACH-00', x = sampleID) & priority == 2) 
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
      colnames(dist_1)[c(1,2,3)] <- c('metasample','sampleID','dist')
      
      if ("Tumors" %in% type) { 
        
        dist_2 <- dist_1 %>% 
          arrange(dist) %>% 
          mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
          filter(grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', x = sampleID) & priority == 2)
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
      
      
      dist_top_n <- dist_2 %>% 
        group_by(metasample) %>%                  
        arrange(dist) %>%  
        slice_head(n = k) %>% 
        ungroup()
      
      top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
      
    }
    
    return(dist_top_n)
    

  } else {
    
    ann_query <- ann[ann$lineage %in% query_lineage,]
    combined_mat <- combined_mat[rownames(combined_mat) %in% ann_query$sampleID,]
    
    if(sum(rownames(combined_mat) %in% input_sample) == 0) {
      showNotification("Select an input lineage")
      warning("Select an input lineage")
      return(NULL)
    }
    
    if(all(!grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', x = rownames(combined_mat)))) {
      showNotification("There are no tumor samples for this lineage in this omics")
      warning("There are no tumor samples for this lineage in this omics")
      return(NULL)
    }
    
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
          filter(grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', x = sampleID) & priority == 2) 
      }
      
      if ("Cell lines" %in% type) {
        
        dist_2 <- dist_1 %>% 
          arrange(dist) %>%
          mutate(priority = if_else(sampleID %in% input_sample, 1, 2)) %>% 
          filter(grepl('ACH-00', x = sampleID) & priority == 2) 
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
      colnames(dist_1)[c(1,2,3)] <- c('metasample','sampleID','dist')
      
      if ("Tumors" %in% type) { 
        
        dist_2 <- dist_1 %>% 
          arrange(dist) %>% 
          mutate(priority = if_else(sampleID %in% selected_samples, 1, 2)) %>% 
          filter(grepl('TCGA|TARGET|TH0|TH1|TH2|TH3|THR', x = sampleID) & priority == 2)
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
      
      
      dist_top_n <- dist_2 %>% 
        group_by(metasample) %>%                  
        arrange(dist) %>%  
        slice_head(n = k) %>% 
        ungroup()
      
      top_k_tumors_1 <- unique(as.character(dist_top_n$sampleID))
      
    }
    
    return(dist_top_n)
    
  }
  
}




