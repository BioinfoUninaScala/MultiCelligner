#' 
#' Get the matrix for CCLE samples and TCGA samples
#'
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import GenomicRanges
#' @import tidyverse
#' @import impute
#' @import MoNETA
#' @import magrittr
#' @import dplyr

############################################################################################ Get TCGA matrix

pancancer_meth <- read.delim("~/celligner/celligner_meth_1/celligner_meth/file/GDC-PANCAN_meth450.tsv.gz", header=FALSE)
id_tcga_pancan <- read_table("file/id_tcga_pancan", col_names = FALSE)
name_pancancer <- c("cpg_id", id_tcga_pancan$X1)
colnames(pancancer_meth) <- name_pancancer

pancancer_meth_ann <- left_join(ann_cpg, pancancer_meth, by = "cpg_id", multiple = "all")
pancancer_meth_ann_sum <- pancancer_meth_ann %>%  group_by(gene) %>% summarise_at(.vars = colnames(pancancer_meth_ann[, -c(1,2)]), .funs = mean)
pancancer_meth_ann_sum_narm <- pancancer_meth_ann %>%  group_by(gene) %>% summarise_at(.vars = colnames(pancancer_meth_ann[, -c(1,2)]), .funs = ~mean(., na.rm = TRUE))

mat_pancancer <- pancancer_meth_ann_sum_narm[, -1]
mat_pancancer <- as.matrix(mat_pancancer)
pancancer_meth_norm <- as.data.frame(mat_pancancer)
rownames(pancancer_meth_norm) <- pancancer_meth_ann_sum$gene

id_pancancer <- colnames(pancancer_meth_norm)
id_pancancer <- substr(id_pancancer, 1, nchar(id_pancancer) -1)
colnames(pancancer_meth_norm) <- id_pancancer

pancancer_meth_ann_index <- which(ann$sampleID %in% colnames(pancancer_meth_norm))
pancancer_ann <- ann[pancancer_meth_ann_index, ]
to_remove_pan <- which(!colnames(pancancer_meth_norm) %in% pancancer_ann$sampleID)
pancancer_meth_norm_filtered <- pancancer_meth_norm[, -to_remove_pan]
pancancer_meth_norm_filtered <- t(pancancer_meth_norm_filtered)

############################################################################################ Get CCLE matrix
################################################################################## Get the cpg_id from RRBS reads

CCLE_RRBS_1kb <- read_delim("CCLE_RRBS_TSS_1kb_20180614.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

loc <- data("Locations")
loc <- as.data.frame(Locations)
loc_range <- GRanges(seqnames = loc$chr, ranges = IRanges(start = loc$pos, end = loc$pos), strand = loc$strand)

rrbs_loc <- CCLE_RRBS_1kb %>% separate(TSS_id, into = c("gene", "chr", "start", "stop"), sep = "_") %>% select("gene", "chr", "start", "stop", "strand")

eliminate_1 <- which(rrbs_loc$chr == "00029157")
eliminate_2 <- which(rrbs_loc$chr == "009911")
rrbs_loc <- rrbs_loc[-eliminate_1, ]
rrbs_loc <- rrbs_loc[-eliminate_2, ]

rrbs_loc$chr <- paste0("chr", rrbs_loc$chr)
rrbs_loc_range <- as(rrbs_loc, "GRanges")

over <- findOverlaps(loc_range, rrbs_loc_range)
merged <- cbind(loc[queryHits(over),], rrbs_loc[subjectHits(over),])
cpg_id <- rownames(merged)

cpg_id <-  gsub("\\.1", "", cpg_id)
cpg_id <- gsub("\\.2", "", cpg_id)
cpg_id <- gsub("\\.3", "", cpg_id)
cpg_id <- gsub("\\.4", "", cpg_id)

merged_1 <- cbind(cpg_id, merged)
ann_cpg <- merged_1[, c("cpg_id","gene")]

rownames(ann_cpg) <- 1:nrow(ann_cpg) 
ann_cpg <- as_tibble(ann_cpg)

############################################################################################ Get CCLE matrix

CCLE_meth_filtered_1 <- left_join(ann_cpg, CCLE_RRBS_1kb, by = "gene", multiple = "all")
CCLE_meth_filtered <- CCLE_meth_filtered_1[, -c(3:8)]

CCLE_meth_filtered_sum <- CCLE_meth_filtered %>%  group_by(gene) %>% summarise_at(.vars = col_to_sum_CCLE, .funs = ~mean(., na.rm = TRUE))

mat_help_na_CCLE <- as.matrix(CCLE_meth_filtered_sum)
rownames(mat_help_na_CCLE) <- mat_help_na_CCLE[, "gene"]
mat_help_na_CCLE <- subset(mat_help_na_CCLE, select = -gene)

na_remove_CCLE <- apply(mat_help_na_CCLE, 1, function (x) all(is.na(x)))
names(na_remove_CCLE) <- NULL
na_remove_CCLE <- which(isTRUE(na_remove_CCLE))

mat_help_na_CCLE_sorted <- mat_help_na_CCLE[, sort(colnames(mat_help_na_CCLE))]
ann_sorted <- ann[order(ann$sampleID_CCLE_Name),]

CCLE_to_remove <- which(!colnames(mat_help_na_CCLE_sorted) %in% ann_sorted$sampleID) 
mat_help_na_CCLE_sorted <- mat_help_na_CCLE_sorted[,-CCLE_to_remove]
CCLE_to_keep_1 <- which(ann_sorted$sampleID %in% colnames(mat_help_na_CCLE_sorted))
CCLE_sample_id <- ann_sorted$sampleID[CCLE_to_keep_1] 

colnames(mat_help_na_CCLE_sorted) <- CCLE_sample_id
meth_mat_CCLE <- t(mat_help_na_CCLE_sorted)
meth_mat_CCLE_1 <- matrix(as.numeric(meth_mat_CCLE), ncol = ncol(meth_mat_CCLE), nrow = nrow(meth_mat_CCLE))
colnames(meth_mat_CCLE_1) <- colnames(meth_mat_CCLE)
rownames(meth_mat_CCLE_1) <- rownames(meth_mat_CCLE)

index_ann_CCLE <- which(ann$sampleID %in% rownames(meth_mat_CCLE_1))
CCLE_meth_ann <- ann[index_ann_CCLE, ]

######################################################################### feature selection

pancancer_meth_norm_filtered_2 <- pancancer_meth_norm_filtered[rowSums(is.na(pancancer_meth_norm_filtered)) < (0.3*ncol(pancancer_meth_norm_filtered)),]
pancancer_meth_norm_filtered_3 <- pancancer_meth_norm_filtered_2[,colSums(is.na(pancancer_meth_norm_filtered)) < (0.8*nrow(pancancer_meth_norm_filtered))]

meth_mat_CCLE_2_1 <- meth_mat_CCLE_1[rowSums(is.na(meth_mat_CCLE_1)) < (0.1*ncol(meth_mat_CCLE_1)),]
meth_mat_CCLE_2_2 <- meth_mat_CCLE_2_1[,colSums(is.na(meth_mat_CCLE_2_1)) < (0.1*nrow(meth_mat_CCLE_2_1))]

lin_remove_1 <- which(pancancer_ann$lineage == "adrenal")
lin_remove_2 <- which(pancancer_ann$lineage == "cervix")
lin_remove_5 <- which(pancancer_ann$lineage == "germ_cell")
lin_remove_6 <- which(pancancer_ann$lineage == "thymus")
tcga_to_remove <- pancancer_ann$sampleID[c(lin_remove_1, lin_remove_2, lin_remove_5, lin_remove_6)]
tcga_to_remove_1 <- which(rownames(pancancer_meth_norm_filtered_3) %in% tcga_to_remove)
pancancer_meth_norm_filtered_3_1 <- pancancer_meth_norm_filtered_3[-tcga_to_remove_1,]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "adrenal"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "cervix"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "germ_cell"),]
pancancer_ann <- pancancer_ann[-which(pancancer_ann$lineage == "thymus"),]

lin_remove_3 <- which(CCLE_meth_ann$lineage == "bone")
lin_remove_4 <- which(CCLE_meth_ann$lineage == "peripheral_nervous_system")
ccle_to_remove <- CCLE_meth_ann$sampleID[c(lin_remove_3, lin_remove_4)]
ccle_to_remove <- which(rownames(meth_mat_CCLE_2_2) %in% ccle_to_remove)
meth_mat_CCLE_2_3 <- meth_mat_CCLE_2_2[-ccle_to_remove,]

CCLE_meth_ann$lineage[which(CCLE_meth_ann$lineage == "plasma_cell")] <- "lymphocyte"
CCLE_meth_ann$lineage[which(CCLE_meth_ann$lineage == "fibroblast")] <- "soft_tissue"
CCLE_meth_ann$lineage[which(CCLE_meth_ann$lineage == "soft tissue")] <- "soft_tissue"
CCLE_meth_ann <- CCLE_meth_ann[-which(CCLE_meth_ann$lineage == "bone"), ]
CCLE_meth_ann <- CCLE_meth_ann[-which(CCLE_meth_ann$lineage == "peripheral_nervous_system"), ]

CCLE_meth_ann_1 <- CCLE_meth_ann[-which(!CCLE_meth_ann$sampleID %in% rownames(meth_mat_CCLE_2_3)),]
meth_mat_CCLE_2_3 <- meth_mat_CCLE_2_3[-which(!rownames(meth_mat_CCLE_2_3) %in% CCLE_meth_ann_1$sampleID),]

################################################################ MATRIX IMPUTATION 

imputed_meth_mat <- NULL

for (j in unique(CCLE_meth_ann_1$lineage)) {
  lineage_vec <- which(CCLE_meth_ann_1$lineage == j)
  lineage_vec_1 <- CCLE_meth_ann_1$sampleID[lineage_vec]
  
  lineage_CCLE_sub_mat <- meth_mat_CCLE_2_3[lineage_vec_1,]
  lineage_CCLE_sub_mat_impute <- impute.knn(lineage_CCLE_sub_mat)
  
  lineage_to_remove <- which(rownames(meth_mat_CCLE_2_3) %in% lineage_vec_1)
  meth_mat_CCLE_2_3 <- meth_mat_CCLE_2_3[-lineage_to_remove, ]
  
  imputed_meth_mat <- rbind(imputed_meth_mat, lineage_CCLE_sub_mat_impute$data)
  
}

CCLE_meth_impute <- rbind(meth_mat_CCLE_2_3, imputed_meth_mat)

################################################################ MATRIX IMPUTATION 

pancancer_meth_norm_filtered_3_1 <- pancancer_meth_norm_filtered_3[-tcga_to_remove_1,]

imputed_meth_mat_tcga <- NULL

for (j in unique(pancancer_ann$lineage)[!unique(pancancer_ann$lineage) %in% c("bile_duct", "nerve")]) {
  lineage_vec_tcga <- which(pancancer_ann$lineage == j)
  lineage_vec_1_tcga <- pancancer_ann$sampleID[lineage_vec_tcga]
  
  lineage_TCGA_sub_mat <- pancancer_meth_norm_filtered_3_1[lineage_vec_1_tcga,]
  
  lineage_TCGA_sub_mat_impute <- impute.knn(lineage_TCGA_sub_mat)
  
  lineage_to_remove_tcga <- which(rownames(pancancer_meth_norm_filtered_3_1) %in% lineage_vec_1_tcga)
  pancancer_meth_norm_filtered_3_1 <- pancancer_meth_norm_filtered_3_1[-lineage_to_remove_tcga, ]
  
  imputed_meth_mat_tcga <- rbind(imputed_meth_mat_tcga, lineage_TCGA_sub_mat_impute$data)
  
}

TCGA_meth_impute <- rbind(pancancer_meth_norm_filtered_3_1, imputed_meth_mat_tcga)

###################################################################### Obtained both matrix:
############################################################## TCGA_meth_impute
############################################################## CCLE_meth_impute


