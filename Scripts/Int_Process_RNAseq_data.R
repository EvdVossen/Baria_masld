rm(list=ls())

# Import Functions & Packages
source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

# imports
ps_rnaseq <- readRDS(file = "raw_data/BARIA.RNAseq.1205.metaV6.220605.RDS")
met <- foreign::read.spss(file = "raw_data/spss epigenetica DS310121.sav", to.data.frame = T) %>% 
  dplyr::mutate(histology = ifelse(histology =="normal", "Normal",
                                   ifelse(histology == "steatosis", "Steatosis", 
                                          "Ballooning")))
subj <- rio::import("Intermediate_files/Subjects_histology.csv")

# gather info
samdat <- data.frame(ps_rnaseq@sam_data) 
samdat_epi <- samdat %>% 
  dplyr::filter(Tissue == "Liver" & Subject_ID %in% subj$Subject_ID) %>%
  dplyr::mutate(Histology = subj$histology[match(.$Subject_ID, subj$Subject_ID)]) %>% 
  dplyr::mutate(Histology = factor(Histology, levels = c("Normal", "Steatosis", "Ballooning"))) %>% 
  `rownames<-`(.$Subject_ID)

# get rnaseq table; take samples with epigenetics info; filter rows (transcripts) with zeros
rnadat <- data.frame(ps_rnaseq@otu_table)
# tax table
taxtab <- data.frame(ps_rnaseq@tax_table)

# First do the normalization step -> Then filter samples that I need.
all(rownames(samdat) == colnames(rnadat))
dds <- DESeqDataSetFromMatrix(countData = round(rnadat),
                              colData = samdat,
                              design = ~ 1)
dds

#Filtering
smallestGroupSize <- 12
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Normalizing factors
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

rnaseq_norm <- as.data.frame(normalized_counts) %>% 
  dplyr::select(all_of(samdat_epi$Sample_id)) %>%
  dplyr::filter(rowSums(.) != 0) %>%
  rename_with(~ samdat_epi$Subject_ID[match(., samdat_epi$Sample_id)], .cols = samdat_epi$Sample_id) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Subject_ID") %>% 
  .[order(as.numeric(.$Subject_ID)),]

rio::export(rnaseq_norm, "Intermediate_files/RNAseq_Normalized_transcripts.csv")