#' 
#' Get multiomics MoNETA similarity matrix
#'
#' @import MoNETA
#' @import dplyr


############################################################################## 

combined_mat <- readRDS("expression_combined_mat.rds") # TCGA and CCLE MNN-corrected expression data matrix
combined_mat_meth <- readRDS("methylation_combined_mat.rds") # TCGA and CCLE MNN-corrected methylation data matrix
combined_mat_mut <- readRDS("mutational_signature_combined_mat.rds") # TCGA and CCLE MNN-corrected mutational signature data matrix

comb_ann_mut <- readRDS("comb_ann_mut.rds") # Mutational signature annotation file
comb_ann_meth <- readRDS("comb_ann_meth.rds") # Methylation annotation file
comb_ann <- readRDS("comb_ann.rds") # Expression annotation file

comb_ann_1 <- comb_ann %>% select(sampleID, lineage, subtype, type)
comb_ann_meth_1 <- comb_ann_meth %>% select(sampleID, lineage, subtype, type)
ann_multiomics <- rbind(comb_ann_mut, comb_ann_1, comb_ann_meth_1)
ann_multiomics <- ann_multiomics %>% select(sampleID, lineage, subtype, type)

############################################################################## multiomics data integration

#  get the list represents the network for each omics with most similar patient connected

multiomics_celligner_knn_MAX_5 <- list(exp_data = k_star_net(t(combined_mat),
                                                             sparsity = .7, 
                                                             distFun = "Euclidean", 
                                                             cores = 60, 
                                                             knn = 5, 
                                                             MAX_ASSOC = 5),
                                       meth_data = k_star_net(t(combined_mat_meth), 
                                                              sparsity = .7, 
                                                              distFun = "Euclidean", 
                                                              cores = 60, 
                                                              knn = 5, 
                                                              MAX_ASSOC = 5),
                                       mut_data = k_star_net(t(combined_mat_mut), 
                                                             sparsity = .7, 
                                                             distFun = "Euclidean", 
                                                             cores = 60, 
                                                             knn = 5, 
                                                             MAX_ASSOC = 5))

# concatenating the individual networks to get multiplex network represents the relationships between 
# nodes across multiple omics

multiplex_knn_MAX_5 <- create_multiplex(multiomics_celligner_knn_MAX_5) 

# start Random Walk with Restart to get patient similarity matrix derived by computing the probability
# distribution of reaching other nodes in the multi-omics network

RWR_mat_knn_MAX_5 <- gen_sim_mat_M(network = multiplex_knn_MAX_5,
                                   tau = NA, restart = 0.7,
                                   jump_neighborhood = F, weighted_multiplex = F, cores = 60)

# dimensionality reduction

emb_RWR_mat_knn_MAX_5_size_32 <- MoNETA::get_embedding(RWR_mat_knn_MAX_5, embedding_size = 32, cores = 50)

emb_umap_RWR_mat_knn_MAX_5_size_32 <- get_parallel_umap_embedding(emb_RWR_mat_knn_MAX_5_size_32, 
                                                                  embedding_size = 2, 
                                                                  n_threads = 60)

emb_tsne_RWR_mat_knn_MAX_5 <- get_tsne_embedding(RWR_mat_knn_MAX_5, 
                                                 embedding_size = 2, 
                                                 perplexity = 70, 
                                                 max_iter = 20000 , 
                                                 num_threads = 80)






