#' 
#' Tuning MNN parameter using Proportion Agreement as alignment quality score
#'

library(readr)
library(tidyverse)
library(MoNETA)
library(magrittr)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(future)
library(furrr)
library(purrr)
library(stats)
library(S4Vectors)
library(foreach)
library(doParallel)
library(batchelor)

source("Celligner_helpers.R")
source("Celligner_method.R")
source("get_prop_agree_v6.R")
source("get_dist_eu_foreach_parallel.R")
source('signature_cluster_data.R')
source('/CreateSeuObj_signature.R')

comb_ann_mut<- readRDS('comb_ann_mut.rds') # TCGA and CCLE annotation file
CCLE_cor_mut<- readRDS('CCLE_cor_mut.rds') # CCLE matrix with the frist 4 cPCs regressed out
TCGA_cor_mut<- readRDS('TCGA_cor_mut.rds') # TCGA matrix with the frist 4 cPCs regressed out
cosmic_signature<- readRDS('cosmic_signature.rds') # Vector of Cosmic signature 
mut_ann_pancancer<- readRDS('TCGA_ann_mut.rds') # TCGA annotation file
mut_ann_ccle<- readRDS("CCLE_ann_mut") # CCLE annotation file

################################# Grid serach of the best combination of parameter (k1,k2,ndist) for the mutual nearest neighbors alignment  
################################################################################## v1

grid_k1 <- seq(5,35, by = 5)
grid_k2 <- seq(5,100, by = 5)
grid_ndist <- 3

grid_mnn_par <- expand.grid(grid_k1, grid_k2, grid_ndist) 
colnames(grid_mnn_par)[c(1,2,3)] <- c("k1", "k2", "ndist")

grid_mnn_par <- t(grid_mnn_par) %>% as.data.frame() 

name_col  <- c()

for (j in 1:ncol(grid_mnn_par)) {
  name_col[j] <- paste(grid_mnn_par[[j]], collapse = "_")
}

colnames(grid_mnn_par) <- name_col

future::plan(multisession, workers = 15) 
plan()

mut_value_prop_mnn_param_v6 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                           CCLE_cor = CCLE_cor_mut,
                                                           TCGA_cor = TCGA_cor_mut,
                                                           TCGA_ann = mut_ann_pancancer,
                                                           CCLE_ann = mut_ann_ccle,
                                                           comb_ann = comb_ann_mut,
                                                           subset_genes = cosmic_signature)


saveRDS(mut_value_prop_mnn_param_v6, "~/celligner/celligner_somatic_mut/grid_search/value_prop_agree/mut_value_prop_mnn_param_v6.rds")

plan("sequential")

##################################################################### v2

grid_k1 <- seq(40,75, by = 5)
grid_k2 <- seq(5,100, by = 5)
grid_ndist <- 3

grid_mnn_par <- expand.grid(grid_k1, grid_k2, grid_ndist) 
colnames(grid_mnn_par)[c(1,2,3)] <- c("k1", "k2", "ndist")

grid_mnn_par <- t(grid_mnn_par) %>% as.data.frame() 

name_col  <- c()

for (j in 1:ncol(grid_mnn_par)) {
  name_col[j] <- paste(grid_mnn_par[[j]], collapse = "_")
}

colnames(grid_mnn_par) <- name_col

future::plan(multisession, workers = 15) 
plan()

mut_value_prop_mnn_param_v6_2 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                             CCLE_cor = CCLE_cor_mut,
                                                             TCGA_cor = TCGA_cor_mut,
                                                             TCGA_ann = mut_ann_pancancer,
                                                             CCLE_ann = mut_ann_ccle,
                                                             comb_ann = comb_ann_mut,
                                                             subset_genes = cosmic_signature)



saveRDS(mut_value_prop_mnn_param_v6_2, "~/celligner/celligner_somatic_mut/grid_search/value_prop_agree/mut_value_prop_mnn_param_v6_2.rds")

plan("sequential")

##################################################################### v3

grid_k1 <- seq(80,115, by = 5)
grid_k2 <- seq(5,100, by = 5)
grid_ndist <- 3

grid_mnn_par <- expand.grid(grid_k1, grid_k2, grid_ndist) 
colnames(grid_mnn_par)[c(1,2,3)] <- c("k1", "k2", "ndist")

grid_mnn_par <- t(grid_mnn_par) %>% as.data.frame() 

name_col  <- c()

for (j in 1:ncol(grid_mnn_par)) {
  name_col[j] <- paste(grid_mnn_par[[j]], collapse = "_")
}

colnames(grid_mnn_par) <- name_col

future::plan(multisession, workers = 15) 
plan()

mut_value_prop_mnn_param_v6_3 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                             CCLE_cor = CCLE_cor_mut,
                                                             TCGA_cor = TCGA_cor_mut,
                                                             TCGA_ann = mut_ann_pancancer,
                                                             CCLE_ann = mut_ann_ccle,
                                                             comb_ann = comb_ann_mut,
                                                             subset_genes = cosmic_signature)



saveRDS(mut_value_prop_mnn_param_v6_3, "~/celligner/celligner_somatic_mut/grid_search/value_prop_agree/mut_value_prop_mnn_param_v6_3.rds")


plan("sequential")

##################################################################### v3

grid_k1 <- seq(120,150, by = 5)
grid_k2 <- seq(5,100, by = 5)
grid_ndist <- 3

grid_mnn_par <- expand.grid(grid_k1, grid_k2, grid_ndist) 
colnames(grid_mnn_par)[c(1,2,3)] <- c("k1", "k2", "ndist")

grid_mnn_par <- t(grid_mnn_par) %>% as.data.frame() 

name_col  <- c()

for (j in 1:ncol(grid_mnn_par)) {
  name_col[j] <- paste(grid_mnn_par[[j]], collapse = "_")
}

colnames(grid_mnn_par) <- name_col

future::plan(multisession, workers = 15) 
plan()

mut_value_prop_mnn_param_v6_4 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                             CCLE_cor = CCLE_cor_mut,
                                                             TCGA_cor = TCGA_cor_mut,
                                                             TCGA_ann = mut_ann_pancancer,
                                                             CCLE_ann = mut_ann_ccle,
                                                             comb_ann = comb_ann_mut,
                                                             subset_genes = cosmic_signature)



saveRDS(mut_value_prop_mnn_param_v6_4, "~/celligner/celligner_somatic_mut/grid_search/value_prop_agree/mut_value_prop_mnn_param_v6_4.rds")


plan("sequential")


############################################################################ close multisession