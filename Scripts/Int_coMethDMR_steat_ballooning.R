rm(list=ls())

#Import packages & functions
source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

# Import data
dat <- get_data(met = T, bval = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)

# #import data
bval <- bval %>% `colnames<-`(gsub("_chr.*","", colnames(.)))

#align metadata to b-value matrix and the comparison of two groups
met_align_steat_ballooning <- as.data.frame(rownames(bval)) %>%
  setNames("Sample") %>% 
  dplyr::mutate(Subject_ID = gsub("BARIA_","", 
                                  gsub("_her OK","", 
                                       gsub("_OK","", 
                                            gsub("replicate_","", 
                                                 gsub("_liver","", Sample)))))) %>% 
  dplyr::mutate(histology = met$histology[match(.$Subject_ID, met$ID)],
                age = met$Age[match(.$Subject_ID, met$ID)],
                sex = met$sex[match(.$Subject_ID, met$ID)],
                bmi = met$bmi[match(.$Subject_ID, met$ID)]) %>% 
  dplyr::mutate(her_ok = ifelse(grepl("her", .$Sample), T, F)) %>% 
  dplyr::filter(histology %in% c("Steatosis", "Ballooning")) %>% 
  dplyr::mutate(age_sq = age ^ 2) %>% 
  dplyr::mutate(histology = as.character(histology),
                is_female = ifelse(sex == 2, TRUE, FALSE),
                is_case = ifelse(histology == "Steatosis", 0, 1)) %>% 
  dplyr::select(Sample, is_female, is_case, age, age_sq, histology) %>% 
  dplyr::mutate(histology = factor(histology, levels = c("Steatosis", "Ballooning")))

bval_steat_ballooning <- bval %>% dplyr::filter(rownames(.) %in% met_align_steat_ballooning$Sample) %>% t(.)

# load pre-defined regions (step 2 of https://github.com/TransBioInfoLab/coMethDMR)
epic_10b4_gene3_200 <- readRDS("Intermediate_files/DMR/EPIC_10b4_Gene_3_200.rds")
epic_10b4_intergene3_200 <- readRDS("Intermediate_files/DMR/EPIC_10b4_InterGene_3_200.rds")

coMethRegpath_steat_ballooning <- "Intermediate_files/DMR/steatosis_ballooning/"
dir.create(coMethRegpath_steat_ballooning, recursive = T)

# residuals -> to remove covariate effects from methylation values by fitting probe-specific linear models and taking the residuals
# resid_obj <- GetResiduals(bval_steat_ballooning,
# betaToM = T,
# pheno_df = met_align_steat_ballooning,
# covariates_char = c("age", "sex", "age_sq"))
# saveRDS(resid_obj, paste0(coMethRegpath_steat_ballooning, "CoMethDMR_files/Residuals_object.RDS"))
resid_obj <- readRDS(paste0(coMethRegpath_steat_ballooning, "CoMethDMR_files/Residuals_object.RDS"))

# Comethallregions -> with adjustment for age sex bmi
# coMethReg <- CoMethAllRegions(dnam = resid_obj,
#                               betaToM = F,
#                               method = "spearman",
#                               rDropThresh_num = 0.4, #recommended according to the article
#                               minCpGs = 3,
#                               genome = "hg19",
#                               arrayType = "EPIC",
#                               CpGs_ls = epic_10b4_gene3_200,
#                               file = NULL,
#                               returnAllCpGs = F,
#                               output = "CpGs",
#                               nCores_int = 10)
# saveRDS(coMethReg, paste0(coMethRegpath_steat_ballooning, "CoMethDMR_files/cometh_residobj_steat_ballooning_gene.RDS"))
coMethReg <- readRDS(paste0(coMethRegpath_steat_ballooning, "CoMethDMR_files/cometh_residobj_steat_ballooning_gene.RDS"))


#intergene
# coMethReg_ig <- CoMethAllRegions(dnam = resid_obj,
#                               betaToM = F,
#                               method = "spearman",
#                               rDropThresh_num = 0.4, #recommended according to the article
#                               minCpGs = 3,
#                               genome = "hg19",
#                               arrayType = "EPIC",
#                               CpGs_ls = epic_10b4_intergene3_200,
#                               file = NULL,
#                               returnAllCpGs = F,
#                               output = "CpGs",
#                               nCores_int = 10)
# saveRDS(coMethReg_ig, paste0(coMethRegpath_steat_ballooning, "cometh_residobj_steat_ballooning_intergene.RDS"))
coMethReg_ig <- readRDS(paste0(coMethRegpath_steat_ballooning, "CoMethDMR_files/cometh_residobj_steat_ballooning_intergene.RDS"))

#for the gene part
# system.time(
#   out_df <- lmmTestAllRegions(
#     betas = bval_steat_ballooning,
#     region_ls = coMethReg,
#     pheno_df = met_align_steat_ballooning,
#     contPheno_char = "is_case",
#     covariates_char = c("is_female", "age", "age_sq"),
#     modelType = "randCoef",
#     arrayType = "EPIC", #change "450k_" to "EPIC_" for EPIC array data
#     nCores_int = 16
#   )
# )
# rio::export(out_df, paste0(coMethRegpath_steat_ballooning, "Result_coMethreg_gene_steat_ballooning.xlsx"))
out_df <- rio::import(paste0(coMethRegpath_steat_ballooning, "Result_coMethreg_gene_steat_ballooning.xlsx"))

#for the intergene part
# system.time(
  # out_df_ig <- lmmTestAllRegions(
#     betas = bval_steat_ballooning,
#     region_ls = coMethReg_ig,
#     pheno_df = met_align_steat_ballooning,
#     contPheno_char = "is_case",
#     covariates_char = c("is_female", "age", "age_sq"),
#     modelType = "randCoef",
#     arrayType = "EPIC", #change "450k_" to "EPIC_" for EPIC array data
#     nCores_int = 16
#   )
# )
# rio::export(out_df_ig, paste0(coMethRegpath_steat_ballooning, "Result_coMethreg_intergene_steat_ballooning.xlsx"))
out_df_ig <- rio::import(paste0(coMethRegpath_steat_ballooning, "Result_coMethreg_intergene_steat_ballooning.xlsx"))

## fix comethreg to dataframe with comma separated cpg's. then make df
cmr_df <- rbind(data.frame(name = names(coMethReg),
                           CpGs = sapply(coMethReg, paste, collapse = ",")),
                data.frame(name = names(coMethReg_ig),
                           CpGs = sapply(coMethReg_ig, paste, collapse = ",")))

# ERROR -> no DMR has a p-value <=0.05.
df <- rbind(out_df %>% filter(FDR <=0.05) %>% mutate(Region = "Gene"),
                out_df_ig %>% filter(FDR <=0.05) %>% mutate(Region = "Intergene")) %>% 
  AnnotateResults(arrayType = "EPIC") %>% 
  mutate(range = paste0(chrom, ":", start, "-", end)) %>%
  dplyr::mutate(CpGs = cmr_df$CpGs[match(range, cmr_df$name)]) %>% 
  .[order(as.numeric(.$FDR)),]

rio::export(df, paste0(coMethRegpath_steat_ballooning, "Result_coMethDMR_Fdr005_steatosis_ballooning.xlsx"))