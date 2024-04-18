#Load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(minfi)
library(coMethDMR)

#Set working directory
setwd("~/Data_files/Epigenetics_Daniela/Data_transfer/")

get_data <- function(met = F, bval = F,  ml_sep_res = F, dmr_sep_res = F, ml_res = F, dmr_res = F, rnaseq_norm = F){
  subj <- rio::import("Intermediate_files/Subjects_histology.csv") %>% 
    dplyr::mutate(histology = factor(histology, levels = c("Normal", "Steatosis", "Ballooning")))
  
  if(met == T){
    met <- foreign::read.spss(file = "raw_data/spss epigenetica DS310121.sav", to.data.frame = T) %>% 
      dplyr::mutate(histology = ifelse(histology =="normal", "Normal",
                                       ifelse(histology == "steatosis", "Steatosis", 
                                              "Ballooning")))
  }
  if(bval == T){
    bval <- rio::import("Intermediate_files/df_beta_val_noob_combat.csv") %>% 
      tibble::column_to_rownames("V1") %>% 
      dplyr::filter(rownames(.) %in% subj$Sample)
  }
  if(ml_sep_res == T){
    ml_res_norm_ball <- rio::import("Intermediate_files/ML/normal_ballooning/output_data/feat_imp_perm_cor_features.xlsx") 
    ml_res_norm_steat <- rio::import("Intermediate_files/ML/normal_steatosis/output_data/feat_imp_perm_cor_features.xlsx") 
    ml_res_steat_ball <- rio::import("Intermediate_files/ML/steatosis_ballooning/output_data/feat_imp_perm_cor_features.xlsx") 
  } else {
    ml_res_norm_ball = F
    ml_res_norm_steat = F
    ml_res_steat_ball = F
  }
  if(dmr_sep_res == T){
    dmr_res_norm_ball <- rio::import("Intermediate_files/DMR/normal_ballooning/Result_coMethDMR_Fdr005_normal_ballooning.xlsx")
    dmr_res_steat_ball <- rio::import("Intermediate_files/DMR/steatosis_ballooning/Result_coMethDMR_Fdr005_steatosis_ballooning.xlsx")
    dmr_three_groups <- rio::import("Intermediate_files/DMR/3_groups/Result_coMethDMR_Fdr005_3groups.xlsx")
  } else{
    dmr_res_norm_ball = F
    dmr_res_steat_ball = F
    dmr_three_groups = F
  }
  # if(ml_res == T){
  #   ml_res <- rio::import("")
  # }
  # if(dmr_res ==T){
  #   dmr_res <- rio::import("")
  # }
  
  if(rnaseq_norm == T){
    ps_rnaseq <- readRDS(file = "raw_data/BARIA.RNAseq.1205.metaV6.220605.RDS")
    taxtab <- data.frame(ps_rnaseq@tax_table)
    rnaseq_norm <- rio::import("Intermediate_files/RNAseq_Normalized_transcripts.csv") %>% 
      tibble::column_to_rownames("Subject_ID") %>% 
      t() %>% 
      as.data.frame()
  } else {
    taxtab = F
  }
  
  d <- list(met, subj, bval, 
            ml_res_norm_ball, ml_res_norm_steat, ml_res_steat_ball, 
            dmr_res_norm_ball, dmr_res_steat_ball, dmr_three_groups,
            rnaseq_norm, taxtab)
  names(d) <- c("met", "subj", "bval", 
                "ml_res_norm_ball", "ml_res_norm_steat", "ml_res_steat_ball", 
                "dmr_res_norm_ball", "dmr_res_steat_ball", "dmr_three_groups",
                "rnaseq_norm", "taxtab")
  
  d <-  base::Filter(function(x) !is.logical(x) || x, d)
  
  # suppressWarnings(dir.create(paste0(path_data, "Manuscript/Main_Figures"), recursive = T))
  # suppressWarnings(dir.create(paste0(path_data, "Manuscript/Supplementary_information"), recursive = T))
  
  return(d)
}
theme_Publication <- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.8), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(0.8)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks.x = element_line(),
            axis.ticks.y = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing  = unit(0, "cm"),
            plot.margin=unit(c(5,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
            plot.caption = element_text(face = "italic", size=rel(0.6))
    ))
} 