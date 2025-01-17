#' 
#' Create a datatable of tumor k nearest neighbors for the metasample sorted by distance, with samples information: lineage and subtype
#'
#' @import dplyr
#' @importFrom SNFtool dist2
#' @importFrom reshape2 melt
#' @param combined_mat combined_mat matrix samples x genes of corrected data by MNN
#' @param reduced_mat dimensionally reduced matrix (tSNE and UMAP): sample x features
#' @param selected_samples vector of samples that will compone the metasample
#' @param n number of the nearest neighbors
#' @param ann annotation file of tumors and cell lines 
#' @return a datatable of tumor samples with lineage and subtype information
#' @export

get_metasample_ann_both_tumor <- function(combined_mat, reduced_mat, selected_samples, n, ann) {
  
x_1 <- combined_mat[which(rownames(combined_mat) %in% selected_samples),]
x_2 <- as.data.frame(colMeans(x_1))
x_3 <- t(x_2)
rownames(x_3) <- 'metasample'
mat_metasample <- as.matrix(x_3)

dist_metasample <- SNFtool::dist2(mat_metasample, combined_mat)
dist_metasmaple_1 <- reshape2::melt(dist_metasample)
colnames(dist_metasmaple_1)[c(1,2,3)] <- c('metasample','sampleID','dist')

dist_top_filtered <- dist_metasmaple_1 %>% filter(grepl('TCGA', sampleID))

dist_top_n <- dist_top_filtered %>%
  group_by(metasample) %>%
  arrange(dist) %>%
  slice_head(n = n + length(selected_samples))

dist_top_n <- dist_top_n[!dist_top_n$sampleID %in% selected_samples,]

dist_top_n_1 <- left_join(dist_top_n[,2], ann[,c(1,2,3)])
dist_top_n_1 <- dist_top_n_1 %>% mutate('n' = 1:nrow(dist_top_n_1))
colnames(dist_top_n_1)[1] <- 'top_neighbors'
dist_top_n_1 <- dist_top_n_1 %>% select(n,top_neighbors,lineage,subtype)
dist_top_n_1 <- dist_top_n_1[!dist_top_n_1$top_neighbors %in% selected_samples,]

return(dist_top_n_1)

}