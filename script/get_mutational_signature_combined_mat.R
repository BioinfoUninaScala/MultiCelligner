#' 
#' Get mutational signature matrix, get mutational signature combined_mat, get_mutational signature confusion matrix and proportion of agreement
#'

library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)
library(readr)
library(tidyverse)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
library(foreach)
library(doParallel)
 

source('signature_cluster_data.R')
source('CreateSeuObj_signature.R')
source('Celligner_helpers.R')
source('Celligner_method.R')
source('Create_Seurat_Object_upload.R')
source('get_heatmap_signature.R')


mut_pancancer_maf <- sigminer::read_xena_variants(path = "~/celligner/celligner_somatic_mut/data/mc3.v0.2.8.PUBLIC.xena.gz")

pancancer_tally <- sig_tally(mut_pancancer_maf,
                             ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                             use_syn = TRUE)

mut_ccle_maf <- sigminer::read_maf(maf = mut_ccle_maf_type_1,
                                   verbose = TRUE)

ccle_tally <- sig_tally(mut_ccle_maf,
                        ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                        use_syn = TRUE)

get_fitting_signature <- function(nmf_mat, cores) {
  
  nmf_mat <- t(nmf_mat)
  x_sign_all <- NULL
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  x_sign_list <- foreach(j = 1:ncol(nmf_mat), .combine = cbind, .packages = "sigminer") %dopar% {
    x_sign <- sig_fit(as.matrix(nmf_mat[,j]), sig_index = "ALL", type = "relative")
    colnames(x_sign) <- colnames(nmf_mat)[j]
    return(x_sign)
  }
  
  stopCluster(cl)
  
  return(x_sign_list)
} # yespar!!!!

mut_fit_sign_ccle <- get_fitting_signature(ccle_tally$nmf_matrix, cores = 30)
mut_fit_sign_pancancer <- get_fitting_signature(pancancer_tally$nmf_matrix, cores = 50)

mut_fit_sign_pancancer_v1 <- mut_fit_sign_pancancer[,which(colnames(mut_fit_sign_pancancer) %in% mut_ann_pancancer$sampleID)]
mut_fit_sign_ccle_v1 <- mut_fit_sign_ccle[,which(colnames(mut_fit_sign_ccle) %in% mut_ann_ccle$ref_ID)]

mut_fit_sign_pancancer_v1 <- t(mut_fit_sign_pancancer_v1)
mut_fit_sign_ccle_v1 <- t(mut_fit_sign_ccle_v1)

TCGA_obj_mut <- createSeuObj_signature(mut_fit_sign_pancancer_v1, mut_ann_pancancer)
CCLE_obj_mut <- createSeuObj_signature(mut_fit_sign_ccle_v1, mut_ann_ccle, type = "CL")

TCGA_obj_mut <- signature_cluster_data(TCGA_obj_mut)
CCLE_obj_mut <- signature_cluster_data(CCLE_obj_mut)

cov_diff_eig_mut <- run_cPCA(TCGA_obj_mut, CCLE_obj_mut, celligner_global$fast_cPCA)

cur_vecs_mut <- cov_diff_eig_mut$rotation[, celligner_global$remove_cPCA_dims, drop = FALSE]
rownames(cur_vecs_mut) <- colnames(mut_fit_sign_pancancer_v1)

TCGA_cor_mut <- resid(lm(t(mut_fit_sign_pancancer_v1) ~ 0 + cur_vecs_mut)) %>% t()
CCLE_cor_mut <- resid(lm(t(mut_fit_sign_ccle_v1) ~ 0 + cur_vecs_mut)) %>% t()

cosmic_signature <- colnames(mut_fit_sign_pancancer_v1)

mnn_res_mut <- run_MNN(CCLE_cor_mut, TCGA_cor_mut,
                       k1 = 20, k2 = 25, ndist = 3,
                       subset_genes = cosmic_signature) 

combined_mat_mut <- rbind(mnn_res_mut$corrected, CCLE_cor_mut)
colnames(mut_ann_ccle)[1] <- "sampleID"
comb_ann_mut <- rbind(mut_ann_pancancer, mut_ann_ccle)
colnames(mut_ann_ccle)[1] <- "ref_ID"
comb_ann_mut <- comb_ann_mut[which(comb_ann_mut$sampleID %in% rownames(combined_mat_mut)),]

comb_obj_mut <- createSeuObj_signature(combined_mat_mut, comb_ann_mut)
comb_obj_mut <- signature_cluster_data(comb_obj_mut)

mnn_param_mut <- c(20,25,3)
prop_agree_mutation <- get_heatmap_signature(mnn_param = mnn_param_mut,
                                            CCLE_cor = CCLE_cor_mut,
                                            TCGA_cor = TCGA_cor_mut,
                                            TCGA_ann = mut_ann_pancancer,
                                            CCLE_ann = mut_ann_ccle,
                                            comb_ann = comb_ann_mut,
                                            subset_genes = cosmic_signature)
