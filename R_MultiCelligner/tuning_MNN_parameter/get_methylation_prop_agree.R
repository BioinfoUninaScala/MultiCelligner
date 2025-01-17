#' 
#' Tuning MNN parameter using Proportion Agreement as alignment quality score
#'
#' @import readr
#' @import tidyverse
#' @import MoNETA
#' @import magrittr
#' @import dplyr
#' @import Seurat
#' @import SeuratObject
#' @import future
#' @import furrr
#' @import purrr
#' @import stats
#' @import S4Vectors
#' @import foreach
#' @import doParallel
#' @import batchelor

source("~/celligner/fun_celligner/celligner_base/Celligner_helpers.R")
source("~/celligner/fun_celligner/celligner_base/Celligner_method.R")
source("~/celligner/fun_celligner/get_prop_agree_v6.R")
source("~/celligner/fun_celligner/get_dist_eu_foreach_parallel.R")

CCLE_cor_meth <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/CCLE_cor_meth.rds")
ccle_meth_ann_1 <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/ccle_meth_ann_1.rds")
comb_ann_meth <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/comb_ann_meth.rds")
DE_gene_set_meth <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/DE_gene_set_meth.rds")
pancancer_ann <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/pancancer_ann.rds")
TCGA_cor_meth <- readRDS("~/celligner/celligner_multiomics/definitive_version/metilazione/script_meth_mat/background_prop_agree/file/TCGA_cor_meth.rds")

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
