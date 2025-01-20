#' 
#' Function to create a confusion matrix (lineage TCGA x lineage CCLE)
#' 
#' @import tidyverse
#' @import dplyr
#' @param combined_mat matrix with MNN corrected data: samples x genes
#' @param ann_multiomics annotation file of all the omics sample
#' @return a list which contains prop_agree scores, confusion matrixs and a dataframe with the top 25 neighbors
#' @export

source("R_MultiCelligner/MultiCelligner_function/get_dist_eu_foreach_parallel.R")

get_confusion_matrix <- function(combined_mat, ann_multiomics) {

dist_1 <- get_fastdist_eu(combined_mat)
print("Dist_couple: did it!")

dist_df <- lapply(dist_1, function(x) do.call(rbind, lapply(x, unlist)))
dist_df <- do.call(rbind, dist_df)
dist_df <- as.data.frame(dist_df)

dist_df_top25 <- dist_df %>% # 
  group_by(CCLEsample) %>%
  arrange(dist_eu) %>%  
  slice_head(n = 25)

colnames(dist_df_top25)[1] <- "ref_ID"
colnames(dist_df_top25)[2] <- "sampleID"

dist_top25_lin <- dplyr::left_join(dist_df_top25, ann_multiomics[, c(1,2)], by = "sampleID") 
colnames(dist_top25_lin)[4] <- "lineage_tcga"

ann_multiomics_1 <- ann_multiomics[grepl("ACH-00", ann_multiomics$sampleID),]
colnames(ann_multiomics_1)[1] <- 'ref_ID'

dist_top25_lin <- dplyr::left_join(dist_top25_lin, ann_multiomics_1[, c(1,2)], by = "ref_ID", ) 
colnames(dist_top25_lin)[5] <- "lineage_ccle"

dist_top25_lin <- dist_top25_lin[,-3]

dist_top25_lin_v2 <- dist_top25_lin %>% 
  group_by(ref_ID) %>%                 
  select(lineage_ccle, lineage_tcga) %>%
  table() %>% as.data.frame()


score_lineage <- dist_top25_lin_v2 %>% group_by(ref_ID) %>% 
  filter(Freq == max(Freq)) %>%                           
  ungroup() %>%
  select(lineage_ccle, lineage_tcga) %>%   

m2_dist <- score_lineage/rowSums(score_lineage) 

m3_dist <- as.data.frame(m2_dist) %>% filter(! is.nan(Freq)) %>% 
  mutate(lineage_ccle = as.character(lineage_ccle), 
         lineage_tcga = as.character(lineage_tcga))

diag_m3_dist <- m3_dist[which(m3_dist$lineage_ccle == m3_dist$lineage_tcga),] 
prop_agree_dist <- sum(diag_m3_dist$Freq) /sum(m3_dist$Freq)               

########################################################

mat_m3_dist <- pivot_wider(m3_dist, names_from = lineage_ccle, values_from = Freq)
mat_m3_dist <- as.matrix(mat_m3_dist)

rownames(mat_m3_dist) <- mat_m3_dist[,1]
mat_m3_dist <- mat_m3_dist[,-1]

mat_m3_dist_1 <- matrix(as.numeric(mat_m3_dist), ncol = ncol(mat_m3_dist), nrow = nrow(mat_m3_dist))

colnames(mat_m3_dist_1) <- colnames(mat_m3_dist)
rownames(mat_m3_dist_1) <- rownames(mat_m3_dist)

########################################################

weig_count <- numeric(nrow(score_lineage))

for (j in 1:nrow(score_lineage)) {
  weig_count[j] <- sum(score_lineage[j,])
}

weig_count <- weig_count / sum(weig_count)

###############################################################

mat_m3_dist_1 <- mat_m3_dist_1 * weig_count

mat_m3_dist_1_d <- as.data.frame(mat_m3_dist_1)

mat_m3_dist_1_d <- mat_m3_dist_1_d %>% mutate(lineage_tcga = rownames(mat_m3_dist_1))

piv_col <- 1:ncol(mat_m3_dist_1_d)

m4_dist <- mat_m3_dist_1_d %>% pivot_longer(cols = piv_col[-length(piv_col)], values_to = "Freq", names_to = "lineage_ccle")

x_dist <- m4_dist[which(m4_dist$lineage_ccle == m4_dist$lineage_tcga),]
prop_agree_weigh_dist <- sum(x_dist$Freq /sum(m4_dist$Freq))

#################################################################################

ccle_heatmap_m3_dist <- ggplot2::ggplot(m3_dist, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
  geom_tile(color = "black") +
  coord_fixed() +
  theme_dark() +
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "CCLE", y = "TCGA", title = "prop_agree_dist") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 

ccle_heatmap_m4_dist <- ggplot2::ggplot(m4_dist, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
  geom_tile(color = "black") +
  coord_fixed() +
  theme_dark() +
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "CCLE", y = "TCGA", title = "prop_agree_weigh_dist") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 

return(list(c(prop_agree_dist, prop_agree_weigh_dist),
            ccle_heatmap_dist = ccle_heatmap_m3_dist,
            ccle_heatmap_weighted_dist = ccle_heatmap_m4_dist,
            dist_df_top25 = dist_df_top25))

}