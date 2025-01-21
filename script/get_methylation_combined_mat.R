#' 
#' Get the methylation combined_mat from the TCGA_meth_impute & CCLE_meth_impute
#'
#' @import Seurat
#' @import SeuratObject
#' @import MoNETA
#' @import magrittr
#' @import dplyr
#' 

source('R_MultiCelligner/celligner_based_function/Celligner_helpers.R')
source('R_MultiCelligner/celligner_based_function/Celligner_method.R')
source('R_MultiCelligner/celligner_based_function/Create_Seurat_Object_upload.R')

common_genes_meth <- colnames(TCGA_meth_impute)

######## remove the lineages which is not present in methylation data
########## remove it form the annotation file and from the data matrix

ccle_meth_ann <- CCLE_meth_ann_1[-which(CCLE_meth_ann_1$lineage == "bile_duct"),]
ccle_meth_ann <- ccle_meth_ann[-which(ccle_meth_ann$lineage == "soft_tissue"),]

pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "soft_tissue"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "eye"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "bile_duct"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "nerve"),]

TCGA_meth_impute_1 <- TCGA_meth_impute[-which(!rownames(TCGA_meth_impute) %in% pancancer_ann$sampleID),]
CCLE_meth_impute_1 <- CCLE_meth_impute[-which(!rownames(CCLE_meth_impute) %in% ccle_meth_ann$sampleID),]

# create seurat objects given an methylation matrix and annotation table

TCGA_obj_meth <- create_Seurat_object(TCGA_meth_impute_1, pancancer_ann, type='tumor')
CCLE_obj_meth <- create_Seurat_object_meth(CCLE_meth_impute_1, ccle_meth_ann, type='CL') 

#take in a Seurat object and run default Seurat clustering algorithm

TCGA_obj_meth <- cluster_data(TCGA_obj_meth)
CCLE_obj_meth <- cluster_data(CCLE_obj_meth)

#calculate gene average expression and variance

gene_stats_meth <- data.frame(
  Tumor_SD = apply(TCGA_meth_impute_1, 2, sd, na.rm=T),
  CCLE_SD = apply(CCLE_meth_impute_1, 2, sd, na.rm=T),
  Tumor_mean = colMeans(TCGA_meth_impute_1, na.rm=T),
  CCLE_mean = colMeans(CCLE_meth_impute_1, na.rm=T),
  Gene = common_genes_meth,
  stringsAsFactors = F) %>% 
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene

#find genes that are differentially expressed between clusters within the data

tumor_DE_genes_meth <- find_differentially_expressed_genes(TCGA_obj_meth) 
CL_DE_genes_meth <- find_differentially_expressed_genes(CCLE_obj_meth) 

DE_genes_meth <- full_join(tumor_DE_genes_meth, CL_DE_genes_meth, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
  dplyr::left_join(gene_stats_meth, by = 'Gene')

#take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors

DE_gene_set_meth <- DE_genes_meth %>%
  dplyr::filter(best_rank < celligner_global$top_DE_genes_per) %>%
  .[['Gene']]

#run contrastive principal components analysis

cov_diff_eig_meth <- run_cPCA(TCGA_obj_meth, CCLE_obj_meth, global$fast_cPCA)

cur_vecs_meth <- cov_diff_eig_meth$rotation[, global$remove_cPCA_dims, drop = FALSE]
rownames(cur_vecs_meth) <- colnames(TCGA_meth_impute_1)

# regressed out the first four cPCs

TCGA_cor_meth <- resid(lm(t(TCGA_meth_impute_1) ~ 0 + cur_vecs_meth)) %>% t()
CCLE_cor_meth <- resid(lm(t(CCLE_meth_impute_1) ~ 0 + cur_vecs_meth)) %>% t()

#run mutual nearest neighbors batch correction

mnn_res_meth <- run_MNN(CCLE_cor_meth, TCGA_cor_meth,
                        k1 = 20, k2 = 15, ndist = 3,
                        subset_genes = DE_gene_set_meth
)

combined_mat_meth <- rbind(mnn_res_meth$corrected, CCLE_cor_meth)

comb_ann_meth <- rbind(pancancer_ann, ccle_meth_ann)

comb_obj_meth <- create_Seurat_object(combined_mat_meth, comb_ann_meth)
comb_obj_meth <- cluster_data(comb_obj_meth)
