#' 
#' Get annotation file with subtype
#'
 
library(tidyverse)
library(dplyr)

############################################################################# create celligner_res with subtype
################################################################## expression sample subtype:

Celligner_res_expression <- Seurat::Embeddings(comb_obj, reduction = "umap") %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
  tibble::rownames_to_column(var = "sampleID") %>%
  dplyr::left_join(comb_obj@meta.data, by = "sampleID")

celligner_res_exp_subtype <- Celligner_res_expression %>% mutate("subtype_1" = paste(lineage, subtype, " "))
celligner_res_exp_subtype <- celligner_res_exp_subtype %>% mutate("subtype_2" = paste(lineage, subtype, " "))

celligner_res_exp_subtype_01 <- celligner_res_exp_subtype[grepl("TCGA", celligner_res_exp_subtype$sampleID),]

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("breast", unique(celligner_res_exp_subtype_01$subtype_1))]  
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast luminal A  ")] <- "breast lumA"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast luminal B  ")] <- "breast lumB"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast HER2-enriched  ")] <- "breast HER2"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast normal  ")] <- "breast"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast NA  ")] <- "breast"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast basal  ")] <- "breast basal"

unique(celligner_res_exp_subtype$subtype_1)[grepl("kidney", unique(celligner_res_exp_subtype$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney kidney clear cell carcinoma  ")] <- "kidney ccRCC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney kidney chromophobe  ")] <- "kidney KICH"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney papillary renal cell carcinoma  ")] <- "kidney PRCC"

unique(celligner_res_exp_subtype$subtype_1)[grepl("skin", unique(celligner_res_exp_subtype$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin neural crest-like  ")] <- "skin NCL"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin melanocytic  ")] <- "skin MEL"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin transitory  ")] <- "skin trans"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin undifferentiated  ")] <- "skin undiff"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin NA  ")] <- "skin"


unique(celligner_res_exp_subtype_01$subtype_1)[grepl("central_nervous_system", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "central_nervous_system glioblastoma multiforme  ")] <- "central_nervous_system GBM"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "central_nervous_system glioma  ")] <- "central_nervous_system LGG"


unique(celligner_res_exp_subtype_01$subtype_1)[grepl("lung", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung lung adenocarcinoma  ")] <- "lung LUAD"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung lung squamous cell carcinoma  ")] <- "lung NSCLC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung NA  ")] <- "lung"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("upper_aerodigestive", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "upper_aerodigestive upper aerodigestive squamous  ")] <- "upper_aerodigestive SCC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("pancreas", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "pancreas adenocarcinoma  ")] <- "pancreas PAAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("uterus", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "uterus uterine carcinosarcoma  ")] <- "uterus UCS"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "uterus uterine endometrioid  ")] <- "uterus UCEC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("thyroid", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "thyroid thyroid carcinoma  ")] <- "thyroid THCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("colorectal", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "colorectal colon adenocarcinoma  ")] <- "colon AC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "colorectal rectum adenocarcinoma  ")] <- "rectum AC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("gastric", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "gastric stomach adenocarcinoma  ")] <- "gastric STAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("ovary", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "ovary ovarian serous cystadenocarcinoma  ")] <- "ovary OV"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("prostate", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "prostate prostate adenocarcinoma  ")] <- "prostate PRAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("urinary_tract", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "urinary_tract bladder urothelial carcinoma  ")] <- "urinary_tract BLCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("liver", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "liver hepatocellular carcinoma  ")] <- "liver LIHC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "liver fibrolamellar hepatocellular carcinoma  ")] <- "liver FLC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("blood", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "blood acute myeloid leukemia  ")] <- "blood LAML"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("esophagus", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "esophagus esophageal carcinoma  ")] <- "esophagus ESCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("lymphocyte", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lymphocyte diffuse large B-cell lymphoma  ")] <- "lymphocyte DLBC"

inf_sub_tcga <- celligner_res_exp_subtype[grepl("TCGA", celligner_res_exp_subtype$sampleID),]
unique(inf_sub_tcga$subtype_2)

celligner_res_exp_subtype_01 <- celligner_res_exp_subtype[grepl("TCGA", celligner_res_exp_subtype$sampleID),]

celligner_res_exp_subtype <- Celligner_res_expression %>% mutate("subtype_1" = paste(lineage, subtype, " "))
celligner_res_exp_subtype <- celligner_res_exp_subtype %>% mutate("subtype_2" = paste(lineage, subtype, " "))

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("breast", unique(celligner_res_exp_subtype_01$subtype_1))]  
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast luminal A  ")] <- "lumA"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast luminal B  ")] <- "lumB"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast HER2-enriched  ")] <- "HER2enr"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast normal  ")] <- "breast norm"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast NA  ")] <- "breast"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "breast basal  ")] <- "breast basal"

unique(celligner_res_exp_subtype$subtype_1)[grepl("kidney", unique(celligner_res_exp_subtype$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney kidney clear cell carcinoma  ")] <- "ccRCC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney kidney chromophobe  ")] <- "KICH"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "kidney papillary renal cell carcinoma  ")] <- "PRCC"

unique(celligner_res_exp_subtype$subtype_1)[grepl("skin", unique(celligner_res_exp_subtype$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin neural crest-like  ")] <- "NCL"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin melanocytic  ")] <- "SKCM"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin transitory  ")] <- "skin trans"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin undifferentiated  ")] <- "skin undiff"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "skin NA  ")] <- "skin"


unique(celligner_res_exp_subtype_01$subtype_1)[grepl("central_nervous_system", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "central_nervous_system glioblastoma multiforme  ")] <- "GBM"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "central_nervous_system glioma  ")] <- "LGG"


unique(celligner_res_exp_subtype_01$subtype_1)[grepl("lung", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung lung adenocarcinoma  ")] <- "LUAD"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung lung squamous cell carcinoma  ")] <- "NSCLC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lung NA  ")] <- "lung"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("upper_aerodigestive", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "upper_aerodigestive upper aerodigestive squamous  ")] <- "SCC"

unique(celligner_res_exp_subtype_01$lineage)

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("pancreas", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "pancreas adenocarcinoma  ")] <- "PAAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("uterus", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "uterus uterine carcinosarcoma  ")] <- "UCS"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "uterus uterine endometrioid  ")] <- "UCEC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("thyroid", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "thyroid thyroid carcinoma  ")] <- "THCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("colorectal", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "colorectal colon adenocarcinoma  ")] <- "colon AC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "colorectal rectum adenocarcinoma  ")] <- "rectum AC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("gastric", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "gastric stomach adenocarcinoma  ")] <- "STAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("ovary", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "ovary ovarian serous cystadenocarcinoma  ")] <- "OV"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("prostate", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "prostate prostate adenocarcinoma  ")] <- "PRAD"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("urinary_tract", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "urinary_tract bladder urothelial carcinoma  ")] <- "BLCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("liver", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "liver hepatocellular carcinoma  ")] <- "LIHC"
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "liver fibrolamellar hepatocellular carcinoma  ")] <- "FLC"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("blood", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "blood acute myeloid leukemia  ")] <- "LAML"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("esophagus", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "esophagus esophageal carcinoma  ")] <- "ESCA"

unique(celligner_res_exp_subtype_01$subtype_1)[grepl("lymphocyte", unique(celligner_res_exp_subtype_01$subtype_1))]
celligner_res_exp_subtype$subtype_2[which(celligner_res_exp_subtype$subtype_2 == "lymphocyte diffuse large B-cell lymphoma  ")] <- "DLBC"

################################################################## methylation sample subtype:

Celligner_res <- Seurat::Embeddings(comb_obj_meth, reduction = "umap") %>%
  as.data.frame() %>%
  magrittr::set_colnames(c("UMAP_1", "UMAP_2")) %>%
  tibble::rownames_to_column(var = "sampleID") %>%
  dplyr::left_join(comb_obj_meth@meta.data, by = "sampleID")

celligner_res_meth_subtype <- Celligner_res_meth %>% mutate("subtype_1" = paste(lineage, subtype, " "))
celligner_res_meth_subtype <- celligner_res_meth_subtype %>% mutate("subtype_2" = paste(lineage, subtype, " "))

celligner_res_meth_subtype_01 <- celligner_res_meth_subtype[grepl("TCGA", celligner_res_meth_subtype$sampleID),]

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("breast", unique(celligner_res_meth_subtype_01$subtype_1))]  
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast luminal A  ")] <- "breast lumA"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast luminal B  ")] <- "breast lumB"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast HER2-enriched  ")] <- "breast HER2"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast normal  ")] <- "breast"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast NA  ")] <- "breast"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast basal  ")] <- "breast basal"



unique(celligner_res_meth_subtype$subtype_1)[grepl("kidney", unique(celligner_res_meth_subtype$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney kidney clear cell carcinoma  ")] <- "kidney ccRCC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney kidney chromophobe  ")] <- "kidney KICH"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney papillary renal cell carcinoma  ")] <- "kidney PRCC"

unique(celligner_res_meth_subtype$subtype_1)[grepl("skin", unique(celligner_res_meth_subtype$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin neural crest-like  ")] <- "skin NCL"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin melanocytic  ")] <- "skin MEL"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin transitory  ")] <- "skin trans"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin undifferentiated  ")] <- "skin undiff"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin NA  ")] <- "skin"


unique(celligner_res_meth_subtype_01$subtype_1)[grepl("central_nervous_system", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "central_nervous_system glioblastoma multiforme  ")] <- "central_nervous_system GBM"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "central_nervous_system glioma  ")] <- "central_nervous_system LGG"


unique(celligner_res_meth_subtype_01$subtype_1)[grepl("lung", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung lung adenocarcinoma  ")] <- "lung LUAD"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung lung squamous cell carcinoma  ")] <- "lung NSCLC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung NA  ")] <- "lung"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("upper_aerodigestive", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "upper_aerodigestive upper aerodigestive squamous  ")] <- "upper_aerodigestive SCC"

unique(celligner_res_meth_subtype_01$lineage)

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("pancreas", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "pancreas adenocarcinoma  ")] <- "pancreas PAAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("uterus", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "uterus uterine carcinosarcoma  ")] <- "uterus UCS"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "uterus uterine endometrioid  ")] <- "uterus UCEC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("thyroid", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "thyroid thyroid carcinoma  ")] <- "thyroid THCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("colorectal", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "colorectal colon adenocarcinoma  ")] <- "colon AC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "colorectal rectum adenocarcinoma  ")] <- "rectum AC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("gastric", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "gastric stomach adenocarcinoma  ")] <- "gastric STAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("ovary", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "ovary ovarian serous cystadenocarcinoma  ")] <- "ovary OV"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("prostate", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "prostate prostate adenocarcinoma  ")] <- "prostate PRAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("urinary_tract", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "urinary_tract bladder urothelial carcinoma  ")] <- "urinary_tract BLCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("liver", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "liver hepatocellular carcinoma  ")] <- "liver LIHC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "liver fibrolamellar hepatocellular carcinoma  ")] <- "liver FLC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("blood", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "blood acute myeloid leukemia  ")] <- "blood LAML"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("esophagus", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "esophagus esophageal carcinoma  ")] <- "esophagus ESCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("lymphocyte", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lymphocyte diffuse large B-cell lymphoma  ")] <- "lymphocyte DLBC"

celligner_res_meth_subtype_01 <- celligner_res_meth_subtype[grepl("TCGA", celligner_res_meth_subtype$sampleID),]

celligner_res_meth_subtype <- Celligner_res_meth %>% mutate("subtype_1" = paste(lineage, subtype, " "))
celligner_res_meth_subtype <- celligner_res_meth_subtype %>% mutate("subtype_2" = paste(lineage, subtype, " "))

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("breast", unique(celligner_res_meth_subtype_01$subtype_1))]  
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast luminal A  ")] <- "lumA"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast luminal B  ")] <- "lumB"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast HER2-enriched  ")] <- "HER2enr"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast normal  ")] <- "breast norm"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast NA  ")] <- "breast"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "breast basal  ")] <- "breast basal"

unique(celligner_res_meth_subtype$subtype_1)[grepl("kidney", unique(celligner_res_meth_subtype$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney kidney clear cell carcinoma  ")] <- "ccRCC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney kidney chromophobe  ")] <- "KICH"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "kidney papillary renal cell carcinoma  ")] <- "PRCC"

unique(celligner_res_meth_subtype$subtype_1)[grepl("skin", unique(celligner_res_meth_subtype$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin neural crest-like  ")] <- "NCL"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin melanocytic  ")] <- "SKCM"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin transitory  ")] <- "skin trans"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin undifferentiated  ")] <- "skin undiff"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "skin NA  ")] <- "skin"


unique(celligner_res_meth_subtype_01$subtype_1)[grepl("central_nervous_system", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "central_nervous_system glioblastoma multiforme  ")] <- "GBM"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "central_nervous_system glioma  ")] <- "LGG"


unique(celligner_res_meth_subtype_01$subtype_1)[grepl("lung", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung lung adenocarcinoma  ")] <- "LUAD"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung lung squamous cell carcinoma  ")] <- "NSCLC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lung NA  ")] <- "lung"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("upper_aerodigestive", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "upper_aerodigestive upper aerodigestive squamous  ")] <- "SCC"

unique(celligner_res_meth_subtype_01$lineage)

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("pancreas", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "pancreas adenocarcinoma  ")] <- "PAAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("uterus", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "uterus uterine carcinosarcoma  ")] <- "UCS"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "uterus uterine endometrioid  ")] <- "UCEC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("thyroid", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "thyroid thyroid carcinoma  ")] <- "THCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("colorectal", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "colorectal colon adenocarcinoma  ")] <- "colon AC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "colorectal rectum adenocarcinoma  ")] <- "rectum AC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("gastric", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "gastric stomach adenocarcinoma  ")] <- "STAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("ovary", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "ovary ovarian serous cystadenocarcinoma  ")] <- "OV"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("prostate", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "prostate prostate adenocarcinoma  ")] <- "PRAD"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("urinary_tract", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "urinary_tract bladder urothelial carcinoma  ")] <- "BLCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("liver", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "liver hepatocellular carcinoma  ")] <- "LIHC"
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "liver fibrolamellar hepatocellular carcinoma  ")] <- "FLC"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("blood", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "blood acute myeloid leukemia  ")] <- "LAML"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("esophagus", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "esophagus esophageal carcinoma  ")] <- "ESCA"

unique(celligner_res_meth_subtype_01$subtype_1)[grepl("lymphocyte", unique(celligner_res_meth_subtype_01$subtype_1))]
celligner_res_meth_subtype$subtype_2[which(celligner_res_meth_subtype$subtype_2 == "lymphocyte diffuse large B-cell lymphoma  ")] <- "DLBC"

