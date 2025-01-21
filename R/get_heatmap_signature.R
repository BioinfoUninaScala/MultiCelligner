#' 
#' Function to generate confusion matrix, Celligner_based_plot, proportion of agreement based on correlation and euclidean distance for mutational signature data
#' 
#' @import tidyverse
#' @import dplyr
#' @import doParallel
#' @import foreach
#' @import readr
#' @import magrittr
#' @import Seurat
#' @import SeuratObject
#' @import stats
#' @import S4Vectors
#' @param mnn_param a vector of value for MNN param (k1,k2,ndist)
#' @param CCLE_cor CCLE matrix with the frist 4 cPCs regressed out
#' @param TCGA_cor TCGA matrix with the frist 4 cPCs regressed out
#' @param TCGA_ann TCGA annotation file
#' @param CCLE_ann CCLE annotation file
#' @param comb_ann TCGA e CCLE annotation file
#' @param subset_genes Vector of genes that have higher variance in both dataset
#' @return a list that contains the confusion matrix, celligner_based plot, proportion of agreement based on correlation and euclidean distance
#' @export

#source("R_MultiCelligner/MultiCelligner_function/get_dist_eu_foreach_parallel.R")
#source("R_MultiCelligner/celligner_based_function/Celligner_method.R")
#source("R_MultiCelligner/celligner_based_function/CreateSeuObj_signature.R")
#source("R_MultiCelligner/celligner_based_function/signature_cluster_data.R")


get_heatmap_signature <- function(mnn_param, CCLE_cor, TCGA_cor, TCGA_ann, CCLE_ann, comb_ann, subset_genes) {
  
  mnn_res <- run_MNN(CCLE_cor, TCGA_cor,
                     k1 = mnn_param[1], 
                     k2 = mnn_param[2], 
                     ndist = mnn_param[3],
                     subset_genes
  )
  
  combined_mat <- rbind(mnn_res$corrected, CCLE_cor)
  
  comb_obj <- createSeuObj_signature(combined_mat, comb_ann)
  comb_obj <- signature_cluster_data(comb_obj)
  
  Celligner_res <- Seurat::Embeddings(comb_obj, reduction = "umap") %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
    tibble::rownames_to_column(var = "sampleID") %>%
    dplyr::left_join(comb_obj@meta.data, by = "sampleID")
  
  lineage_averages <- Celligner_res %>%
    dplyr::filter(!lineage %in% c(
      "embryo", "endocrine", "engineered", "engineered_blood",
      "engineered_breast", "engineered_central_nervous_system", "engineered_kidney",
      "engineered_lung", "engineered_ovary", "engineered_prostate", "epidermoid_carcinoma",
      "nasopharynx", "nerve", "pineal", "teratoma", "unknown"
    )) %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1, na.rm = T),
      UMAP_2 = median(UMAP_2, na.rm = T)
    )
  lineage_averages$lineage <- gsub("_", " ", lineage_averages$lineage)
  lineage_lab_aes <- ggplot2::geom_text(data = lineage_averages, mapping = aes(x = UMAP_1, y = UMAP_2, label = lineage), size = 3, color = "#000000")
  
  celligner_plot <- ggplot2::ggplot(Celligner_res, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(alpha = 0.7, pch = 21, ggplot2::aes(color = type, fill = lineage, size = type)) +
    ggplot2::scale_color_manual(values = c(tumor = "white", CL = "black")) +
    ggplot2::scale_size_manual(values = c(tumor = 0.75, CL = 1.5)) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::guides(
      fill = "none",
      color = ggplot2::guide_legend(override.aes = list(color = c("black", "white"), fill = c("white", "black")))
    ) +
    ggplot2::theme_classic()
  
  print("UMAP: did it!")
  
  tumor_CL_cor <- calc_tumor_CL_cor(combined_mat, comb_ann)
  
  print("tumor_CL_cor: did it!")
  
  oneccle_topcor_25tcga <- list()
  
  for (i in 1:ncol(tumor_CL_cor)) {
    top_25_rows <- head(order(tumor_CL_cor[, i], decreasing = TRUE), 25)
    
    col_choose <- tumor_CL_cor[top_25_rows, i]
    
    df_ccle_cor <- data.frame(col_choose)
    
    nam_colonna <- colnames(tumor_CL_cor)
    colonna_orig_nam <- nam_colonna[i]
    colnames(df_ccle_cor) <- colonna_orig_nam
    
    colnames(df_ccle_cor) <- "pearson_correlation"
    
    df_ccle_cor$sampleID <- rownames(df_ccle_cor)
    
    df_ccle_cor <- mutate(ref_ID = rep(colonna_orig_nam, 25), df_ccle_cor)
    
    df_ccle_cor <- dplyr::left_join(df_ccle_cor, TCGA_ann[, c(1,3)], by = "sampleID")
    colnames(df_ccle_cor)[4] <- "lineage_tcga"
    
    df_ccle_cor <- dplyr::left_join(df_ccle_cor, CCLE_ann[, c(1,3)], by = "ref_ID")
    colnames(df_ccle_cor)[5] <- "lineage_ccle"
    
    oneccle_topcor_25tcga[[colonna_orig_nam]] <- df_ccle_cor
  }
  
  table_oneccle_topcor_25tcga <- oneccle_topcor_25tcga
  
  lineage_col <- "lineage_tcga"
  
  table_oneccle_topcor_25tcga <- list()
  
  for (nome_df in names(oneccle_topcor_25tcga)) {
    
    df <- oneccle_topcor_25tcga[[nome_df]]
    
    conteggio <- table(df[[lineage_col]])
    
    as_df <- as.data.frame(conteggio)
    
    as_df$Freq <- as_df$Freq / 25
    
    table_oneccle_topcor_25tcga[[nome_df]] <- as_df
    
  }
  
  ######################################
  
  table_oneccle_topcor_25tcga_ccle <- oneccle_topcor_25tcga
  
  lineage_col <- "lineage_ccle"
  
  table_oneccle_topcor_25tcga_ccle <- list()
  
  for (nome_df in names(oneccle_topcor_25tcga)) {
    
    df <- oneccle_topcor_25tcga[[nome_df]]
    
    conteggio <- table(df[[lineage_col]])
    
    as_df <- as.data.frame(conteggio)
    
    as_df$Freq <- as_df$Freq / 25 
    
    table_oneccle_topcor_25tcga_ccle[[nome_df]] <- as_df
    
  }
  
  mat_freq_heatmap <- bind_rows(table_oneccle_topcor_25tcga, .id = "sampleID")
  
  mat_freq_heatmap_1 <- bind_rows(table_oneccle_topcor_25tcga_ccle, .id = "sampleID")
  
  mat_freq_heatmap_all <- left_join(mat_freq_heatmap, mat_freq_heatmap_1, by = "sampleID" )
  
  colnames(mat_freq_heatmap_all)[c(2,4)] <- c("lineage_tcga", "lineage_ccle")
  colnames(mat_freq_heatmap_all)[c(3,5)] <- c("freq_tcga", "freq_ccle")
  
  m1_cor <- mat_freq_heatmap_all %>% group_by(sampleID) %>% 
    filter(freq_tcga == max(freq_tcga)) %>% 
    ungroup() %>%
    select(lineage_ccle, lineage_tcga) %>% 
    table 
  
  
  m2_cor <- m1_cor/rowSums(m1_cor) 
  
  m3_cor <- as.data.frame(m2_cor) %>% filter(! is.nan(Freq)) %>% 
    mutate(lineage_ccle = as.character(lineage_ccle), 
           lineage_tcga = as.character(lineage_tcga))
  
  diag_m3_cor <- m3_cor[which(m3_cor$lineage_ccle == m3_cor$lineage_tcga),]
  prop_agree_cor <- sum(diag_m3_cor$Freq) /sum(m3_cor$Freq)
  
  ####################################################################
  
  mat_m3_cor <- pivot_wider(m3_cor, names_from = lineage_ccle, values_from = Freq) 
  mat_m3_cor <- as.matrix(mat_m3_cor)
  
  rownames(mat_m3_cor) <- mat_m3_cor[,1]
  mat_m3_cor <- mat_m3_cor[,-1]
  
  mat_m3_cor_1 <- matrix(as.numeric(mat_m3_cor), ncol = ncol(mat_m3_cor), nrow = nrow(mat_m3_cor))
  
  colnames(mat_m3_cor_1) <- colnames(mat_m3_cor)
  rownames(mat_m3_cor_1) <- rownames(mat_m3_cor)
  
  ########################################################
  
  weig_count <- numeric(nrow(m1_cor))
  
  for (j in 1:nrow(m1_cor)) {
    weig_count[j] <- sum(m1_cor[j,])
  }
  
  weig_count <- weig_count / sum(weig_count)
  
  ###############################################################
  
  mat_m3_cor_1 <- mat_m3_cor_1 * weig_count
  
  mat_m3_cor_1_d <- as.data.frame(mat_m3_cor_1)
  
  mat_m3_cor_1_d <- mat_m3_cor_1_d %>% mutate(lineage_ccle = rownames(mat_m3_cor_1))
  
  piv_col <- 1:ncol(mat_m3_cor_1_d)
  
  m4_cor <- mat_m3_cor_1_d %>% pivot_longer(cols = piv_col[-length(piv_col)], values_to = "Freq", names_to = "lineage_tcga")
  
  x <- m4_cor[which(m4_cor$lineage_ccle == m4_cor$lineage_tcga),]
  prop_agree_weigh_cor <- sum(x$Freq /sum(m4_cor$Freq))
  
  ######################################################## finish cor
  
  print("Correlation: did it!")
  
  ######################################################## start dist_eu
  
  dist_couple <- get_fastdist_eu(mnn_res = mnn_res, CCLE_cor = CCLE_cor, cores = 35)
  
  print("Dist_couple: did it!")
  
  dist_df <- lapply(dist_couple, function(x) do.call(rbind, lapply(x, unlist)))
  rm(dist_couple)
  dist_df <- do.call(rbind, dist_df)
  dist_df <- as.data.frame(dist_df)
  
  dist_df_top25 <- dist_df %>%
    group_by(CCLEsample) %>%
    arrange(dist_eu) %>%  
    slice_head(n = 25)
  
  rm(dist_df)
  
  colnames(dist_df_top25)[1] <- "ref_ID"
  colnames(dist_df_top25)[2] <- "sampleID"
  
  dist_top25_lin <- dplyr::left_join(dist_df_top25, TCGA_ann[, c(1,3)], by = "sampleID")
  colnames(dist_top25_lin)[4] <- "lineage_tcga"
  
  dist_top25_lin <- dplyr::left_join(dist_top25_lin, CCLE_ann[, c(1,3)], by = "ref_ID", )
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
    table                                  
  
  rm(dist_top25_lin_v2)
  
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
  
  print("im doing the heatmap")
  
  ccle_heatmap_m3_cor <- ggplot2::ggplot(m3_cor, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
    geom_tile(color = "black") +
    coord_fixed() +
    theme_dark() +
    scale_fill_gradient(low = "white", high = "red") + 
    labs(x = "CCLE", y = "TCGA", title = "prop_agree_cor") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
    guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 
  
  ccle_heatmap_m4_cor <- ggplot2::ggplot(m4_cor, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
    geom_tile(color = "black") +
    coord_fixed() +
    theme_dark() +
    scale_fill_gradient(low = "white", high = "red") + 
    labs(x = "CCLE", y = "TCGA", title = "prop_agree_weigh_cor") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
    guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 
  
  ccle_heatmap_m3_dist <- ggplot2::ggplot(m3_dist, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
    geom_tile(color = "black") +
    coord_fixed() +
    theme_dark() +
    scale_fill_gradient(low = "white", high = "red") + 
    labs(x = "CCLE", y = "TCGA", title = "prop_agree_dist") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
    guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 
  
  ccle_heatmap_m4_dist <- ggplot2::ggplot(m4_dist, aes(x = lineage_ccle, y = lineage_tcga, fill = Freq)) +
    geom_tile(color = "black") +
    coord_fixed() +
    theme_dark() +
    scale_fill_gradient(low = "white", high = "red") + 
    labs(x = "CCLE", y = "TCGA", title = "prop_agree_weigh_dist") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #ruota i nomi sull'asse x
    guides(fill = guide_colourbar(barwidth = 0.5,barheight = 20)) 
  
  print("Heatmap: did it!")
  
  #################################################################################
  
  return(list(c(prop_agree_cor = prop_agree_cor, prop_agree_weigh_cor = prop_agree_weigh_cor, 
                prop_agree_dist = prop_agree_dist, prop_agree_weigh_dist = prop_agree_weigh_dist),
              ccle_heatmap_cor = ccle_heatmap_m3_cor,
              ccle_heatmap_weighted_cor = ccle_heatmap_m4_cor,
              ccle_heatmap_dist = ccle_heatmap_m3_dist,
              ccle_heatmap_weighted_dist = ccle_heatmap_m4_dist,
              celligner_plot = celligner_plot + lineage_lab_aes))
  
}
