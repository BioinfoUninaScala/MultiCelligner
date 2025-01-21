#' 
#' Get expression combined_mat
#'
#' @import Seurat
#' @import SeuratObject
#' @import MoNETA
#' @import magrittr
#' @import dplyr
#' @import readr
#' @import tidyverse
#' 

source('R_MultiCelligner/celligner_based_function/Celligner_helpers.R')
source('R_MultiCelligner/celligner_based_function/Celligner_method.R')
source('R_MultiCelligner/celligner_based_function/Create_Seurat_Object_upload.R')
source('R_MultiCelligner/MultiCelligner_function/get_heatmap_umap_v6.R')

TCGA_mat_exp <- readRDS("~/celligner/celligner_19Q4_exp/data/TCGA_mat_exp.rds")
hgnc_complete_set <- readRDS("~/celligner/celligner_19Q4_exp/data/hgnc_complete_set.rds")

common_genes <- intersect(colnames(TCGA_mat_exp), hgnc_complete_set$Gene)
TCGA_mat_exp <- TCGA_mat_exp[,common_genes]
hgnc_complete_set <- filter(hgnc_complete_set, Gene %in% common_genes)
hgnc_complete_set <- hgnc_complete_set[!duplicated(hgnc_complete_set$Gene),]
rownames(hgnc_complete_set) <- hgnc_complete_set$Gene
colnames(TCGA_mat_exp) <- hgnc_complete_set$Gene

CCLE_mat_1 <- read_csv("~/celligner/celligner_19Q4_exp/data/CCLE_expression_full.csv")

CCLE_mat <- CCLE_mat_1 %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('...1') %>%
  as.matrix()

colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]

Celligner_info <- read.csv("~/celligner/celligner_19Q4_exp/data/Celligner_info.csv", header=FALSE)
colnames(Celligner_info) <- Celligner_info[1, ]
Celligner_info <- Celligner_info[-1, ]

Celligner_info_sub <- Celligner_info[-which(grepl("TARGET", Celligner_info$sampleID)),]
Celligner_info_sub <- Celligner_info_sub[-which(grepl("TH", Celligner_info_sub$sampleID)),]

ann <- Celligner_info_sub %>% as.data.frame()
if('UMAP_1' %in% colnames(ann)) {
  ann <- ann %>% 
    dplyr::select(-UMAP_1)
}
if('UMAP_2' %in% colnames(ann)) {
  ann <- ann %>% 
    dplyr::select(-UMAP_2)
}
if('cluster' %in% colnames(ann)) {
  ann <- ann %>% 
    dplyr::select(-cluster)
}

TCGA_ann <- dplyr::filter(ann, type=='tumor')
CCLE_ann <- dplyr::filter(ann, type=='CL')

######## remove the lineages which is not present in other omics data
########## remove it form the annotation file and from the data matrix

TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "nerve"),] 
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "germ_cell"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "thymus"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "soft_tissue"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "eye"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "cervix"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "adrenal"),]
TCGA_ann <- TCGA_ann[-which(TCGA_ann$lineage == "bile_duct"),]

CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "soft_tissue"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "fibroblast"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "bile_duct"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "cervix"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "bone"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "peripheral_nervous_system"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "eye"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "adrenal_cortex"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_prostate"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_kidney"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_central_nervous_system"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "embryo"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_ovary"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_lung"),] 
CCLE_ann <- CCLE_ann[-which(CCLE_ann$lineage == "engineered_breast"),] 

CCLE_ann$lineage[which(CCLE_ann$lineage == "plasma_cell")] <- "lymphocyte"

x1 <- rownames(TCGA_mat_exp)
x1_1 <- gsub("\\.", "-", x1)
rownames(TCGA_mat_exp) <- x1_1

TCGA_mat_exp <- TCGA_mat_exp[which(rownames(TCGA_mat_exp) %in% TCGA_ann$sampleID),]

CCLE_mat <- CCLE_mat[which(rownames(CCLE_mat) %in% CCLE_ann$sampleID),]

genes_used <- intersect(colnames(TCGA_mat_exp), colnames(CCLE_mat))

TCGA_mat <- TCGA_mat_exp[,genes_used]
CCLE_mat <- CCLE_mat[,genes_used]

colnames(CCLE_mat) <- hgnc_complete_set$Symbol
colnames(TCGA_mat) <- hgnc_complete_set$Symbol

common_genes <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))

gene_stats <- data.frame(
  Tumor_SD = apply(TCGA_mat, 2, sd, na.rm=T),
  CCLE_SD = apply(CCLE_mat, 2, sd, na.rm=T),
  Tumor_mean = colMeans(TCGA_mat, na.rm=T),
  CCLE_mean = colMeans(CCLE_mat, na.rm=T),
  Gene = common_genes,
  stringsAsFactors = F) %>% 
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) 

colnames(hgnc_complete_set)[c(1,2)] <- c("Symbol", "Gene")
gene_stats <- left_join(hgnc_complete_set, gene_stats, by = "Gene")

comb_ann <- rbind(
  TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'tumor'),
  CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'CL')
)

TCGA_obj <- create_Seurat_object(TCGA_mat, TCGA_ann, type='tumor')
CCLE_obj <- create_Seurat_object(CCLE_mat, CCLE_ann, type='CL')

TCGA_obj <- cluster_data(TCGA_obj)
CCLE_obj <- cluster_data(CCLE_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
CL_DE_genes <- find_differentially_expressed_genes(CCLE_obj)

DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
  dplyr::left_join(gene_stats, by = 'Gene')

DE_gene_set <- DE_genes %>%
  dplyr::filter(best_rank < celligner_global$top_DE_genes_per) %>%
  .[['Gene']]

cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, celligner_global$fast_cPCA)

cur_vecs <- cov_diff_eig$rotation[, celligner_global$remove_cPCA_dims, drop = FALSE]
rownames(cur_vecs) <- colnames(TCGA_mat)

TCGA_cor <- resid(lm(t(TCGA_mat) ~ 0 + cur_vecs)) %>% t()
CCLE_cor <- resid(lm(t(CCLE_mat) ~ 0 + cur_vecs)) %>% t()

mnn_res <- run_MNN(CCLE_cor, TCGA_cor,  k1 = celligner_global$mnn_k_tumor, k2 = celligner_global$mnn_k_CL, ndist = celligner_global$mnn_ndist, 
                   subset_genes = DE_gene_set)

combined_mat <- rbind(mnn_res$corrected, CCLE_cor)

comb_obj <- create_Seurat_object(combined_mat, comb_ann)
comb_obj <- cluster_data(comb_obj)

mnn_param_exp <- c(celligner_global$mnn_k_tumor, celligner_global$mnn_k_CL, celligner_global$mnn_ndist)

prop_agree_expression <- get_heatmap_umap_v6(mnn_param = mnn_param_exp, 
                                             CCLE_cor = CCLE_cor,
                                             TCGA_cor = TCGA_cor,
                                             TCGA_ann = TCGA_ann,
                                             CCLE_ann = CCLE_ann_1,
                                             comb_ann = comb_ann,
                                             subset_genes = DE_gene_set)

