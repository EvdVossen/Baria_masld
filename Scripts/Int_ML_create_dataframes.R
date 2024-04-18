rm(list=ls())

#Import packages & functions
source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

# Import data
dat <- get_data(met = T, bval = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)

#align metadata to beta-value matrix
met_align <- as.data.frame(rownames(bval)) %>%
  setNames("Sample") %>% 
  dplyr::mutate(Subject_ID = gsub("BARIA_","", 
                                       gsub("_OK","", 
                                            gsub("_liver","", Sample)))) %>% 
  dplyr::mutate(histology = met$histology[match(.$Subject_ID, met$ID)])

#make 3 separate groups
met_steat <- met_align %>% dplyr::filter(histology == "Steatosis")
met_norm <- met_align %>% dplyr::filter(histology == "Normal")
met_ball <- met_align %>% dplyr::filter(histology =="Ballooning")

##making the groups (taking the 8000 features with the largest variance)
#normal vs ballooning
xtrain_n_b <- bval %>% 
  dplyr::filter(rownames(.) %in% c(met_norm$Sample, met_ball$Sample)) %>% 
  dplyr::select(all_of(names(sort(sapply(., var), decreasing = TRUE)[1:8000]))) %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID, everything())

ytrain_n_b <- xtrain_n_b %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID) %>% 
  dplyr::mutate(CLASS_Good = ifelse(rownames(.) %in% met_ball$Sample, 1,0))

#normal vs steatosis
xtrain_n_s <- bval %>% 
  dplyr::filter(rownames(.) %in% c(met_norm$Sample, met_steat$Sample)) %>% 
  dplyr::select(all_of(names(sort(sapply(., var), decreasing = TRUE)[1:8000]))) %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID, everything())

ytrain_n_s <- xtrain_n_s %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID) %>% 
  dplyr::mutate(CLASS_Good = ifelse(rownames(.) %in% met_steat$Sample, 1,0))

#steatosis vs ballooning
xtrain_s_b <- bval %>% 
  dplyr::filter(rownames(.) %in% c(met_steat$Sample, met_ball$Sample)) %>% 
  dplyr::select(all_of(names(sort(sapply(., var), decreasing = TRUE)[1:8000]))) %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID, everything())

ytrain_s_b <- xtrain_s_b %>% 
  dplyr::mutate(SAMPLE_ID = rownames(.)) %>% 
  dplyr::select(SAMPLE_ID) %>% 
  dplyr::mutate(CLASS_Good = ifelse(rownames(.) %in% met_ball$Sample, 1,0))

##Save dataframes
#create paths
ml_path <- "Intermediate_files/ML/"
dir.create(paste0(ml_path, "normal_ballooning"), recursive = T)
dir.create(paste0(ml_path, "normal_steatosis"))
dir.create(paste0(ml_path, "steatosis_ballooning"))

# #normal vs ballooning
rio::export(xtrain_n_b, paste0(ml_path, "normal_ballooning/xtrain.xlsx"))
rio::export(ytrain_n_b, paste0(ml_path, "normal_ballooning/ytrain.xlsx"))

#normal vs steatosis
rio::export(xtrain_n_s, paste0(ml_path, "normal_steatosis/xtrain.xlsx"))
rio::export(ytrain_n_s, paste0(ml_path, "normal_steatosis/ytrain.xlsx"))

#steatosis vs ballooning
rio::export(xtrain_s_b, paste0(ml_path, "steatosis_ballooning/xtrain.xlsx"))
rio::export(ytrain_s_b, paste0(ml_path, "steatosis_ballooning/ytrain.xlsx"))