rm(list=ls())

#Import packages & functions
source("~/Data_files/Epigenetics_Daniela/Data_transfer/R-scripts/functions.R")

# Import data
dat <- get_data(met = T, bval = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)


#import files
bval <- rio::import("noob_combat_noXY_Snps_GAPsnps/df_beta_val_noob_combat.csv") %>% 
  tibble::column_to_rownames("V1") %>% 
  filter(!grepl("her OK", rownames(.)))
met <- foreign::read.spss(file = "../spss epigenetica DS310121.sav", to.data.frame = T)

#align metadata to beta-value matrix
met_align <- as.data.frame(rownames(bval)) %>%
  setNames("Sample") %>% 
  dplyr::mutate(Subject_ID = gsub("BARIA_","", 
                                  gsub("_her OK","", 
                                       gsub("_OK","", 
                                            gsub("replicate_","", 
                                                 gsub("_liver","", Sample)))))) %>% 
  dplyr::mutate(histology = met$histology[match(.$Subject_ID, met$ID)])

#make 3 separate groups
met_steat <- met_align %>% filter(histology == "steatosis")
met_norm <- met_align %>% filter(histology == "normal")
met_ball <- met_align %>% filter(histology =="balooning")

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
ml_path <- "ML/three_separate_models/"
dir.create(paste0(ml_path, "normal_ballooning"), recursive = T)
dir.create(paste0(ml_path, "normal_steatosis"))
dir.create(paste0(ml_path, "steatosis_ballooning"))

# # #normal vs ballooning
# rio::export(xtrain_n_b, paste0(ml_path, "normal_ballooning/xtrain.xlsx"))
# rio::export(ytrain_n_b, paste0(ml_path, "normal_ballooning/ytrain.xlsx"))
# 
# #normal vs steatosis
# rio::export(xtrain_n_s, paste0(ml_path, "normal_steatosis/xtrain.xlsx"))
# rio::export(ytrain_n_s, paste0(ml_path, "normal_steatosis/ytrain.xlsx"))
# 
# #steatosis vs ballooning
# rio::export(xtrain_s_b, paste0(ml_path, "steatosis_ballooning/xtrain.xlsx"))
# rio::export(ytrain_s_b, paste0(ml_path, "steatosis_ballooning/ytrain.xlsx"))


##Check on collinearity
n_b_cor <- rcorr(as.matrix(xtrain_n_b %>% select(-SAMPLE_ID)))
test <- data.frame("X1" = rownames(n_b_cor$r)[row(n_b_cor$r)], 
                   "X2" = colnames(n_b_cor$r)[col(n_b_cor$r)],
                   "rho" = c(n_b_cor$r),
                   "p_value" = c(n_b_cor$P))
xtrain_n_b1 <- xtrain_n_b %>% 
  dplyr::mutate(CLASS_Good = ytrain_n_b$CLASS_Good[match(xtrain_n_b$SAMPLE_ID, ytrain_n_b$SAMPLE_ID)])
ggscatter(data = xtrain_n_b1, x = "cg00063654", y = "cg27540865",
          color = "CLASS_Good",
          add = "reg.line",               # Add regression line
          conf.int = TRUE,                # Show confidence interval
          cor.coef = TRUE,                # Show correlation coefficient
          cor.method = "pearson",        # Use Pearson correlation coefficient
          # xlab = "cg24760581",
          # ylab = "cg17480035",
          title = "Scatterplot with Regression Line and Confidence Interval") +
  theme_minimal()


test2 <- FSinR::featureSelection(data = xtrain_n_b1, 
                                 class = "CLASS_Good",
                                 searcher = ,
                                 evaluator = "fscore")

