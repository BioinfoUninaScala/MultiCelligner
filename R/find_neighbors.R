#' 
#' Create a data frame with the neighbors and related annotation of a query sample/s
#'
#' @import dplyr
#' @import magrittr
#' @importFrom SNFtool dist2
#' @importFrom reshape2 melt
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param input_sample a vector of one or multiple TCGA and/or CCLE IDs. 
#' @param k number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @param query_type type of the k neighbors samples (Tumor or Cell Line)
#' @param query_lineage limit the neighbors searching to that lineage/s
#' @return an interactive plot that highlighting the tumor k nearest neighbors 
#' @export
#' 

find_neighbors <- function(combined_mat, reduced_mat, input_sample, k, ann, query_type, query_lineage = 'All') {

  if ("All" %in% query_lineage && length(query_lineage) > 1) {
    shiny::showNotification("Select only All or choose multiple lineage without All")
    warning("Select only All or choose multiple lineage without All")
    return(NULL)
  }
  
  ### when you switch omics with the sample, if the there isn't the sample in that omics will appear a notification!
  if(all(!input_sample %in% colnames(reduced_mat))) {
    shiny::showNotification("There isn't this input sample(s) in this omics")
    warning("There isn't this input sample(s) in this omics")
    return(NULL)
  }
  
  ann <- ann %>%
    dplyr::mutate(type = case_when(
      type == "Tumor" ~ "Tumors", 
      type == "CL" ~ "Cell lines")
    ) %>% dplyr::filter(type %in% query_type)
  
  if(all(query_lineage == 'All')) {
    ann_query <- ann
  } else {
    ann_query <- ann[ann$lineage %in% query_lineage,]
  }
  
  ann_query2 <- union(ann_query$sampleID, input_sample)
  
  
  combined_mat <- combined_mat[rownames(combined_mat) %in% ann_query2,]
  
  if(all(!input_sample %in% rownames(combined_mat))) {
    shiny::showNotification("There isn't this input sample/s for this LINEAGE in this omics")
    warning("There isn't this input sample/s for this LINEAGE in this omics")
    return(NULL)
  }
  
  x_1 <- combined_mat[which(rownames(combined_mat) %in% input_sample), , drop = FALSE]
  
  if(length(input_sample) == 1) {
    x_3 <- x_1
  }else{
    x_2 <- as.data.frame(colMeans(x_1))
    x_3 <- t(x_2)
    rownames(x_3) <- 'metasample'
    x_3 <- as.matrix(x_3)
  }
  
  query_combined_mat <- combined_mat[which(!rownames(combined_mat) %in% input_sample),]
  dist <- SNFtool::dist2(x_3, query_combined_mat)
  dist_1 <- reshape2::melt(dist) %>% dplyr::mutate_if(is.factor, as.character)
  colnames(dist_1)[c(1,2,3)] <- c('ref_ID', 'sampleID', 'dist')
  
  dist_top_n <- dist_1 %>% 
    dplyr::arrange(dist)  %>% 
    dplyr::group_by(ref_ID) %>%                  
    dplyr::arrange(dist) %>%  
    dplyr::slice_head(n = k) %>% 
    dplyr::ungroup()
  
  return(dist_top_n)
  
}
