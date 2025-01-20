#' 
#' Get multiomics Similarity Network Fusion (SNF) matrix
#'
#' @import SNFtool
#' @import dplyr
#' @import tidyverse

### Calculate the pair-wise distance;

SNF_dist_exp_2 <- SNFtool::dist2(as.matrix(mofa_mat_exp), as.matrix(mofa_mat_exp)) # expression matrix with omics intersected samples
SNF_dist_meth_2 <- SNFtool::dist2(as.matrix(mofa_mat_meth), as.matrix(mofa_mat_meth)) # methylation matrix with omics intersected samples
SNF_dist_mut_2 <- SNFtool::dist2(as.matrix(mofa_mat_mut), as.matrix(mofa_mat_mut)) # mutational signature matrix with omics intersected samples

### construct similarity graphs

SNF_w1 <- affinityMatrix(SNF_dist_exp_2)
SNF_w2 <- affinityMatrix(SNF_dist_meth_2)
SNF_w3 <- affinityMatrix(SNF_dist_mut_2)

### fuse all the graphs
### then the overall matrix can be computed by similarity network fusion (SNF):

SNF_W <- SNFtool::SNF(list(SNF_w1,SNF_w2,SNF_w3))
colnames(SNF_W) <- colnames(SNF_w1)
rownames(SNF_W) <- rownames(SNF_w1)




