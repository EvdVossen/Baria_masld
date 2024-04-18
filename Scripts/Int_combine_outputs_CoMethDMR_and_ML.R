rm(list=ls())

source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

dat <- get_data(met = T, ml_sep_res = T, dmr_sep_res = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)

#Combine DMR lists
rbind(dmr_res_norm_ball, dmr_res_steat_ball, dmr_res_three_groups) %>% 
  .[order(as.numeric(.$pValue)),] %>% 
  rio::export(., "Intermediate_files/DMR/Sig_DMRs.xlsx")

#Combine top 10 from ML
bind_rows((ml_res_norm_ball %>% 
             head(10) %>%
             dplyr::mutate(AUC = 0.86,
                           ML_run = "Normal v Ballooning",
                           Group_1 = "Normal",
                           Group_2 = "Ballooning")),
          (ml_res_steat_ball %>% 
             head(10) %>% 
             dplyr::mutate(AUC = 0.82,
                           ML_run = "Steatosis v Ballooning",
                           Group_1 = "Steatosis",
                           Group_2 = "Ballooning")),
          (ml_res_norm_steat %>% 
             head(10) %>% 
             dplyr::mutate(AUC = 0.83,
                           ML_run = "Normal v Steatosis",
                           Group_1 = "Normal",
                           Group_2 = "Steatosis"))) %>%
  rio::export(., "Intermediate_files/ML/Top_ML_results.xlsx")
