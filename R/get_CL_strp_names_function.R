#' 
#' Get the stripped cell line name of the user selected combined matrix
#' 
#' @import dplyr
#' @param combined_mat matrix samples x genes of corrected data by MNN: 
#' @param ann annotation file of tumors and cell lines 
#' @return a vector with the stripped cell line name
#' @export


get_CL_strp_names <- function(combined_mat ,ann) {
  
  depmap_id <- rownames(combined_mat)
  depmap_id_df <- as.data.frame(depmap_id)
  colnames(depmap_id_df)[1] <- 'sampleID'
  df <- left_join(depmap_id_df, ann[,c(1,6)], by = 'sampleID')
  nm <- df$stripped_cell_line_name[!grepl('TCGA', df$stripped_cell_line_name)]
  
  return(nm)
  
}