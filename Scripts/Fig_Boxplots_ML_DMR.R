rm(list=ls())

source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

dat <- get_data(met = T, dmr_res = T, ml_res = T, bval = T)
base::list2env(dat,envir=.GlobalEnv); rm(dat)

#Remove additional information of the CpG in the column name
bval1 <- bval %>% `colnames<-`(gsub("_chr.*","", colnames(.)))


#Generate ML results
ml_res_loop <- pivot_longer(ml_res, cols = c(FeatName, correlated_features), names_to = "Variable", values_to = "CpG_full") %>%
  dplyr::mutate(CpG = gsub("_chr.*","", CpG_full)) %>%
  dplyr::filter(!is.na(CpG_full)) %>%
  tidyr::separate_rows(CpG_full, sep = ",") %>%
  dplyr::mutate(CpG = gsub("_chr.*","", CpG_full)) %>%
  dplyr::select(CpG, CpG_full, RelFeatImp, Group_1, Group_2)

#all CpGs that were found in the DMRs
dmr_cpgs <- dmr_res %>% 
  mutate(strings = strsplit(.$CpGs, ",")) %>% 
  tidyr::unnest(strings) %>% 
  .$strings %>% 
  .[!duplicated(.)]


## Reduce the search space for bval matrix (time-efficiency)
bval_ml_dmr <- bval1 %>% select(all_of(c(dmr_cpgs, ml_res_loop$CpG))) %>% 
  dplyr::mutate(Subject_ID = subj$Subject_ID[match(rownames(.), subj$Sample)],
                Histology = subj$histology[match(rownames(.), subj$Sample)])

#Create directory for saving the results
suppressWarnings(dir.create("Plots/Boxplots/ML", recursive = T))

# loop to create the boxplots for ML (separate)
for (idx in 1:nrow(ml_res_loop)){
  #Select the index including the CpG
  CpG_dat <- ml_res_loop[idx,]
  
  #Isolate CpG name for convenience of making the boxplot
  name_cpg <- CpG_dat$CpG
  
  #Isolate the CpG
  df_cpg <- bval_ml_dmr %>% 
    dplyr::select(all_of(name_cpg), Subject_ID, Histology) %>% 
    setNames(c("name_cpg", "Subject_ID", "Histology"))
  
  #Create the boxplot
  p_l <- 
    ggplot(df_cpg, aes(x = Histology, 
                       y = name_cpg, 
                       fill = Histology)) +
    geom_boxplot(alpha = 0.5, 
                 width = 0.4, 
                 outlier.shape = NA) +
    geom_jitter(aes(col = name_cpg), 
                size = 1, 
                color = "black", 
                width = 0.2, 
                height = 0) +
    ylab(label = CpG_dat$CpG_full) +
    xlab(label = "Histology") +
    theme_Publication() +
    scale_fill_manual(name = "Histology", 
                      labels = c("Normal", "Steatosis", "Ballooning"),
                      values = c(pal_hist[1], pal_hist[2], pal_hist[3])) +
    theme(legend.text.align = 0,
          strip.background = element_rect(fill = NA, 
                                          color = NA)) + 
    stat_compare_means(comparisons = list(c("Normal", "Steatosis"), 
                                          c("Normal", "Ballooning"), 
                                          c("Steatosis", "Ballooning")), 
                       method = "wilcox.test",
                       label = "p.format", 
                       paired = F)  +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  #get the rounded feature important to easily trace back in the map the importance
  feat_imp_rounded <- round(CpG_dat$RelFeatImp)
  
  #Save boxplot
  ggsave(filename = paste0(feat_imp_rounded, "_", 
                           CpG_dat$Group_1, "_", 
                           CpG_dat$Group_2, "_", 
                           CpG_dat$CpG_full, ".png"), 
         plot = p_l, 
         device = "png",
         path = "Plots/Boxplots/ML", 
         width = 6, 
         height = 6)
}

# loop to create the boxplots for DMRs (all CpGs per DMR)
for (idx in 1:nrow(dmr_res)){
  #Set filename to save to either the name of the gene the CpG is located; for intergenes put the intergene + its chromosome position with range
  filename = ifelse(!is.na(dmr_res[idx,]$UCSC_RefGene_Name), 
                    dmr_res[idx,]$UCSC_RefGene_Name,
                    paste0(dmr_res[idx,]$Region, "_", dmr_res[idx,]$range))
  
  cat(paste0("loop number: ", idx, "\nname_gene or intergene DMR location: ", filename, "\n\n"))
  
  #Isolate CpGs of the DMR of the loop
  dmr_cpg_oi = dmr_res[idx,] %>% 
    mutate(strings = strsplit(.$CpGs, ",")) %>% 
    tidyr::unnest(strings) %>% 
    .$strings %>% 
    .[!duplicated(.)]
  
  #Select the isolated CpGs
  df_cpg <- bval_ml_dmr %>% 
    dplyr::select(all_of(dmr_cpg_oi), Subject_ID, Histology) 
  
  #New loop to create the boxplots separately
  pl <- list()
  for (idxj in 1:length(dmr_cpg_oi)){
    #Isolate CpG
    cpg_name = dmr_cpg_oi[idxj]
    
    #Create df with one CpG
    df_cpg_fig <- df_cpg %>% 
      dplyr::select(all_of(cpg_name), Subject_ID, Histology) %>% 
      setNames(c("name_cpg", "Subject_ID", "Histology"))
    
    #Create the boxplot
    pl[[idxj]] <- 
      ggplot(df_cpg_fig, aes(x = Histology, 
                             y = name_cpg, 
                             fill = Histology)) +
      geom_boxplot(alpha = 0.5, 
                   width = 0.4, 
                   outlier.shape = NA) +
      geom_jitter(aes(col = name_cpg), 
                  size = 1, 
                  color = "black", 
                  width = 0.2, 
                  height = 0) +
      ylab(label = CpG_dat$CpG_full) +
      xlab(label = "Histology") +
      theme_Publication() +
      scale_fill_manual(name = "Histology", 
                        labels = c("Normal", "Steatosis", "Ballooning"),
                        values = c(pal_hist[1], pal_hist[2], pal_hist[3])) +
      theme(legend.text.align = 0,
            strip.background = element_rect(fill = NA, 
                                            color = NA)) + 
      stat_compare_means(comparisons = list(c("Normal", "Steatosis"), 
                                            c("Normal", "Ballooning"), 
                                            c("Steatosis", "Ballooning")), 
                         method = "wilcox.test",
                         label = "p.format", 
                         paired = F)  +
      guides(fill = guide_legend(override.aes = list(shape = NA)))
    
  }
  
  #wrap boxplots
  p_cpg_dmr <- patchwork::wrap_plots(pl, ncol = ifelse(length(pl) > 5, 3, 2)) + plot_layout(guides = "collect")
  
  #save the CpGs of a DMR together 
  ggsave(filename = paste0(filename,".png"), plot = p_cpg_dmr, device = "png",
         path = "Plots/Boxplots/DMR", 
         width = ifelse(length(pl) > 6, 16, 12),
         height = ifelse(length(pl) > 6, 12, 8))
}
