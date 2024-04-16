rm(list=ls())

# Import Functions & Packages
source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

# Import data
dat <- get_data(met = T, bval = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)

# Create an aligned df with subject information
met_align <- as.data.frame(rownames(bval)) %>%
  setNames("Sample") %>% 
  dplyr::mutate(Subject_ID = gsub("BARIA_","", 
                                  gsub("_her OK","", 
                                       gsub("_OK","", 
                                            gsub("replicate_","", 
                                                 gsub("_liver","", Sample)))))) %>% 
  merge(x = ., y = met, by.x = "Subject_ID", by.y = "ID") %>% 
  dplyr::select(Sample, Subject_ID, histology)
rio::export(met_align, "Intermediate_files/Subjects_histology.csv")
