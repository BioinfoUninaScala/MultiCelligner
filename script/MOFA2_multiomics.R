#' 
#' Get multiomics MOFA2 matrix
#'
#' @import MOFA2
#' @import magrittr
#' @import dplyr
#' @import tidyverse
#' @import gplots
#' 

############################   MOFA INTEGRATION ################################

################# intersect the omics samples

combined_mat <- readRDS('combined_mat.rds')
combined_mat_meth <- readRDS('combined_mat_meth.rds') # Methylation matrix of MNN-corrected data
combined_mat_mut <- readRDS('combined_mat_mut.rds')

omics_sample <- intersect(rownames(combined_mat), rownames(combined_mat_meth))
omics_sample_1 <- intersect(omics_sample, rownames(combined_mat_mut))

mofa_mat_exp <- combined_mat[which(rownames(combined_mat) %in% omics_sample_1),]
mofa_mat_meth <- combined_mat_meth[which(rownames(combined_mat_meth) %in% omics_sample_1),]
mofa_mat_mut <- combined_mat_mut[which(rownames(combined_mat_mut) %in% omics_sample_1),]

mofa_mat_exp <- mofa_mat_exp[omics_sample_1,]
mofa_mat_meth <- mofa_mat_meth[omics_sample_1,]
mofa_mat_mut <- mofa_mat_mut[omics_sample_1,]

multi_celligner <- list(exp = t(mofa_mat_exp),
                        meth = t(mofa_mat_meth),
                        mut = t(mofa_mat_mut))

MOFAobject <- create_mofa(multi_celligner)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"MOFA_multiomics_model.hdf5")

Sys.setenv(OMP_NUM_THREADS=100) ### do it in parallel

MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

model <- load_model('MOFA_multiomics_model.hdf5')

mofa_ann <- ann_multiomics[ann_multiomics$sampleID %in% omics_sample_1, ]
mofa_ann <- mofa_ann[-which(duplicated(mofa_ann$sampleID)),]

sample_metadata <- mofa_ann %>% dplyr::select(sample=sampleID, dplyr::everything())
all(samples_names(MOFAobject)[[1]] == sample_metadata$sample)

samples_metadata(model) <- sample_metadata

model <- run_tsne(model)
model <- run_umap(model)

MOFA <- model@expectations$Z$group1

###  Dimensionality Reduction

dta_mofa_umap <- model@dim_red$UMAP[,-1] %>% t()
all(colnames(dta_mofa_umap) == mofa_ann$sampleID)
MOFA_UMAP <-plot_2D_matrix(dta_mofa_umap, nodes_anno = mofa_ann, id_name = "sampleID", interactive = TRUE,
                           id_anno_color = "lineage", id_anno_shape = "type",
                           wo_legend = FALSE, title = "MOFA UMAP")


dta_mofa_tsne <- model@dim_red$TSNE[,-1] %>% t()
colnames(dta_mofa_tsne) <- model@dim_red$TSNE$sample
all(colnames(dta_mofa_tsne) == mofa_ann$sampleID)
MOFA_TSNE <- plot_2D_matrix(dta_mofa_tsne , nodes_anno = mofa_ann, id_name = "sampleID", interactive = TRUE,
                            id_anno_color = "lineage", id_anno_shape = "type",
                            wo_legend = FALSE, title = "MOFA t-SNE")

