#' 
#' Get the stripped cell line name and the tumors ID of the user selected combined matrix
#' 
#' @import dplyr
#' @param combined_mat matrix samples x genes of corrected data by MNN 
#' @param ann annotation file of tumors and cell lines 
#' @return a vector with the stripped cell line name and the tumors ID
#' @export

get_all_name_strpp <- function(combined_mat, ann, sample_info) {
  
  nm <- rownames(combined_mat)
  nm_1 <- as.data.frame(nm)
  colnames(nm_1) <- 'sampleID'
  
  s_inf <- sample_info[,c(1,3)]
  colnames(s_inf)[1] <- 'sampleID'

  x <- left_join(nm_1, s_inf, by = 'sampleID')
  
  x$stripped_cell_line_name[which(is.na(x$stripped_cell_line_name))] <- x$sampleID[grepl('TCGA', x$sampleID)]
  
  x$stripped_cell_line_name[which(x$sampleID == 'ACH-001173')] <- 'CVCL_VU83'
  x$stripped_cell_line_name[which(x$sampleID == 'ACH-001316')] <- 'ACH-001316'
  x$stripped_cell_line_name[which(x$sampleID == 'ACH-000010')] <- 'ACH-000010'
  
  return(x$stripped_cell_line_name)
  
}