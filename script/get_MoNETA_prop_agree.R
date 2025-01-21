#' 
#' Get Proportion of Agreement score for multiomics MoNETA similarity matrix 
#'
#' @importFrom reshape2 melt
#' @import MoNETA
#' @import magrittr
#' @import dplyr
#' @import tidyverse
#' @import ggplot2
#' 

####################################################################### get_similarity matrix

RWR_mat_knn_MAX_5 <- gen_sim_mat_M(network = multiplex_knn_MAX_5,
                                   tau = NA, restart = 0.7,
                                   jump_neighborhood = F, weighted_multiplex = F, cores = 80)

####################################################################### calculate proportion of agreement 

prop_agree_RWR_df_1 <- reshape2::melt(RWR_mat_knn_MAX_5)
prop_agree_RWR_df_2 <- prop_agree_RWR_df_1[which(grepl("ACH-", prop_agree_RWR_df_1$Var1)),]
prop_agree_RWR_df_2 <- prop_agree_RWR_df_2[-which(grepl("TCGA", prop_agree_RWR_df_1$Var1)),]
prop_agree_RWR_df_3 <- prop_agree_RWR_df_2[-which(grepl("ACH", prop_agree_RWR_df_2$Var2)),]
prop_agree_RWR_df <- prop_agree_RWR_df_3 %>% group_by(Var1) %>%  arrange(desc(value)) %>% slice_head(n = 25)
prop_agree_RWR_df <- prop_agree_RWR_df[-which(grepl("TCGA", prop_agree_RWR_df$Var1)),]


colnames(prop_agree_RWR_df)[1] <- c("sampleID")
prop_RWR_df_v1 <- left_join(prop_agree_RWR_df, ann_multiomics[c(1,2)], by = "sampleID")

colnames(prop_RWR_df_v1)[4] <- "lineage_ccle"

colnames(prop_agree_RWR_df)[c(1,2)] <- c("Var1", "sampleID")

prop_RWR_df_help <- left_join(prop_agree_RWR_df[,2], ann_multiomics[c(1,2)], by = "sampleID")
colnames(prop_RWR_df_help)[c(1,2)] <- c("Var2","lineage_tcga")
prop_RWR_df_help <- prop_RWR_df_help[-which(duplicated(prop_RWR_df_help$Var2)),]

prop_RWR_df_v2 <- left_join(prop_RWR_df_v1, prop_RWR_df_help, by = "Var2")

prop_RWR_df_v2 <- prop_RWR_df_v2[-which(duplicated(prop_RWR_df_v2$value)),]

##############################################################################################

multiplex_score_lineage_1 <- 
  prop_RWR_df_v2 %>% group_by(sampleID) %>% 
  select(lineage_ccle, lineage_tcga) %>%
  table() %>% as.data.frame()

multiplex_score_lineage <- multiplex_score_lineage_1 %>% 
  group_by(sampleID) %>% 
  dplyr::filter(Freq == max(Freq)) %>%
  ungroup() %>% 
  select(lineage_ccle, lineage_tcga) %>% 
  table 

m2_multiplex <- multiplex_score_lineage/rowSums(multiplex_score_lineage) 

m3_multipelx <- as.data.frame(m2_multiplex) %>% dplyr::filter(! is.nan(Freq)) %>% 
  mutate(lineage_ccle = as.character(lineage_ccle), 
         lineage_tcga = as.character(lineage_tcga))

diag_m3_multiplex <- m3_multipelx[which(m3_multipelx$lineage_ccle == m3_multipelx$lineage_tcga),] 
prop_agree_dist <- sum(diag_m3_multiplex$Freq) /sum(m3_multipelx$Freq)

##############################################################################################

ccle_heatmap_m3_multiplex <- ggplot2::ggplot(m3_multipelx, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
  geom_tile(color = "black") +
  coord_fixed() +
  theme_dark() +
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "CCLE", y = "TCGA", title = "prop_agree_dist") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 


##############################################################################################

mat_m3_multiplex <- pivot_wider(m3_multipelx, names_from = lineage_ccle, values_from = Freq)
mat_m3_multiplex <- as.matrix(mat_m3_multiplex)

rownames(mat_m3_multiplex) <- mat_m3_multiplex[,1]
mat_m3_multiplex <- mat_m3_multiplex[,-1]

mat_m3_multiplex_1 <- matrix(as.numeric(mat_m3_multiplex), ncol = ncol(mat_m3_multiplex), nrow = nrow(mat_m3_multiplex))

colnames(mat_m3_multiplex_1) <- colnames(mat_m3_multiplex)
rownames(mat_m3_multiplex_1) <- rownames(mat_m3_multiplex)

##############################################################################################

weig_count <- numeric(nrow(multiplex_score_lineage))

for (j in 1:nrow(multiplex_score_lineage)) {
  weig_count[j] <- sum(multiplex_score_lineage[j,])
}

weig_count <- weig_count / sum(weig_count)

##############################################################################################

mat_m3_multiplex_1 <- mat_m3_multiplex_1 * weig_count

mat_m3_multiplex_1_d <- as.data.frame(mat_m3_multiplex_1)

mat_m3_multiplex_1_d <- mat_m3_multiplex_1_d %>% mutate(lineage_ccle = rownames(mat_m3_multiplex_1))

piv_col <- 1:ncol(mat_m3_multiplex_1_d)

m4_multiplex <- mat_m3_multiplex_1_d %>% pivot_longer(cols = piv_col[-length(piv_col)], values_to = "Freq", names_to = "lineage_tcga")

diag_m4_multiplex <- m4_multiplex[which(m4_multiplex$lineage_ccle == m4_multiplex$lineage_tcga),]
prop_agree_weigh_dist <- sum(diag_m4_multiplex$Freq /sum(m4_multiplex$Freq))

ccle_heatmap_m4_multiplex <- ggplot2::ggplot(m4_multiplex, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
  geom_tile(color = "black") +
  coord_fixed() +
  theme_dark() +
  scale_fill_gradient(low = "white", high = "red") + 
  labs(x = "CCLE", y = "TCGA", title = "prop_agree_weigh_dist") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
  guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 
