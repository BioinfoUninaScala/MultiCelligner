#' 
#' Tuning MNN parameter using Proportion Agreement as alignment quality score
#'

library(readr)
library(tidyverse)
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

CCLE_cor_meth <- readRDS('CCLE_cor_meth.rds') # CCLE matrix with the frist 4 cPCs regressed out
ccle_meth_ann_1 <- readRDS('CCLE_ann_meth.rds') # CCLE annotation file
comb_ann_meth <- readRDS('comb_ann_meth.rds') # TCGA e CCLE annotation file
DE_gene_set_meth <- readRDS('DE_gene_set_meth.rds') # Vector of genes that have higher variance in both dataset
pancancer_ann <- readRDS('TCGA_ann_meth.rds') # TCGA annotation file
TCGA_cor_meth <- readRDS('TCGA_cor_meth.rds') # TCGA matrix with the frist 4 cPCs regressed out

################################# Grid serach of the best combination of parameter (k1,k2,ndist) for the mutual nearest neighbors alignment  
################################################################################## v1

grid_k1 <- seq(5,40, by = 5)
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

future::plan(multisession, workers = 17) 
plan()

value_prop_mnn_param_v6_1 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                         CCLE_cor = CCLE_cor_meth,
                                                         TCGA_cor = TCGA_cor_meth,
                                                         TCGA_ann = pancancer_ann,
                                                         CCLE_ann = ccle_meth_ann_1,
                                                         comb_ann = comb_ann_meth,
                                                         subset_genes = DE_gene_set_meth)


saveRDS(value_prop_mnn_param_v6_1, "~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/result_file/value_prop_mnn_param_v6_1.rds")

plan("sequential")

############################################################### v2


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

future::plan(multisession, workers = 17) 
plan()

value_prop_mnn_param_v6_2 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                         CCLE_cor = CCLE_cor_meth,
                                                         TCGA_cor = TCGA_cor_meth,
                                                         TCGA_ann = pancancer_ann,
                                                         CCLE_ann = ccle_meth_ann_1,
                                                         comb_ann = comb_ann_meth,
                                                         subset_genes = DE_gene_set_meth)


saveRDS(value_prop_mnn_param_v6_2, "~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/result_file/value_prop_mnn_param_v6_2.rds")

plan("sequential")

############################################################### v3

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

future::plan(multisession, workers = 17) 
plan()

value_prop_mnn_param_v6_3 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                         CCLE_cor = CCLE_cor_meth,
                                                         TCGA_cor = TCGA_cor_meth,
                                                         TCGA_ann = pancancer_ann,
                                                         CCLE_ann = ccle_meth_ann_1,
                                                         comb_ann = comb_ann_meth,
                                                         subset_genes = DE_gene_set_meth)


saveRDS(value_prop_mnn_param_v6_3, "~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/result_file/value_prop_mnn_param_v6_3.rds")

plan("sequential")

################################################################## v4

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

future::plan(multisession, workers = 17) 
plan()

value_prop_mnn_param_v6_4 <- grid_mnn_par %>% future_map(get_prop_agree_v6, 
                                                         CCLE_cor = CCLE_cor_meth,
                                                         TCGA_cor = TCGA_cor_meth,
                                                         TCGA_ann = pancancer_ann,
                                                         CCLE_ann = ccle_meth_ann_1,
                                                         comb_ann = comb_ann_meth,
                                                         subset_genes = DE_gene_set_meth)


saveRDS(value_prop_mnn_param_v6_4, "~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/result_file/value_prop_mnn_param_v6_4.rds")


plan("sequential")
