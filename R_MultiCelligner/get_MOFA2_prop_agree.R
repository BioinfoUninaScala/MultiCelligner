#' 
#' Get Proportion of Agreement score for multiomics MOFA2 matrix 
#'
#' @import dplyr
#' @import tidyverse
#' 

source("R_MultiCelligner/MultiCelligner_function/get_MOFA_fastdist_eu.R")

MOFA <- readRDS("MOFA_mat.rds") # load the MOFA2 matrix
ann_multiomics <- readRDS("ann_multiomics.rds") # load multiomics annotation file

dist_couple <- get_MOFA_fastdist_eu(MOFA)

#######################################################################################

dist_df <- lapply(dist_couple, function(x) do.call(rbind, lapply(x, unlist)))

dist_df <- do.call(rbind, dist_df)
dist_df <- as.data.frame(dist_df)

dist_df_top25 <- dist_df %>% 
  group_by(CCLEsample) %>%
  arrange(dist_eu) %>%  
  slice_head(n = 25)

colnames(dist_df_top25)[1] <- "ref_ID"
colnames(dist_df_top25)[2] <- "sampleID"

dist_top25_lin <- dplyr::left_join(dist_df_top25, ann_multiomics[, c(1,2)], by = "sampleID")
colnames(dist_top25_lin)[4] <- "lineage_tcga"

ann_multiomics_1 <- ann_multiomics
colnames(ann_multiomics_1)[1] <- "ref_ID"

dist_top25_lin <- dplyr::left_join(dist_top25_lin, ann_multiomics_1[, c(1,2)], by = "ref_ID", )
colnames(dist_top25_lin)[5] <- "lineage_ccle"


dist_top25_lin <- dist_top25_lin[,-3]
dist_top25_lin <- dist_top25_lin[-which(duplicated(dist_top25_lin$sampleID)),]

dist_top25_lin_v2 <- dist_top25_lin %>% 
  group_by(ref_ID) %>%                
  select(lineage_ccle, lineage_tcga) %>%
  table() %>% as.data.frame()


score_lineage <- dist_top25_lin_v2 %>% group_by(ref_ID) %>% 
  dplyr::filter(Freq == max(Freq)) %>%                        
  ungroup() %>%
  select(lineage_ccle, lineage_tcga) %>% 
  table                                  

m2_dist <- score_lineage/rowSums(score_lineage) 

m3_dist <- as.data.frame(m2_dist) %>% dplyr::filter(! is.nan(Freq)) %>% 
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

