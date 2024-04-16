rm(list=ls())

source("~/Data_files/Epigenetics_Daniela/Data_transfer/Baria_masld/Scripts/functions.R")

library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(limma)
library(affy)
library(pheatmap)
library(wateRmelon)

workdirectory=paste0(getwd(), "/raw_data/liver_dna_methylation/")
setwd(workdirectory)

idatspath=paste0(workdirectory, "idats/") 

files <- list.files(paste0(workdirectory, "STS_files"))

file1 <- utils::read.table(paste0(workdirectory,"STS_files/", files[2]), sep=",", skip = 7, header=T) %>% 
  dplyr::mutate(Batch = "File1")

file2 <- utils::read.table(paste0(workdirectory,"STS_files/", files[3]), sep=",", skip = 8, header=T) %>% 
  dplyr::mutate(Batch = "File2") %>% 
  dplyr::select(names(.)[names(.) %in% names(file1)])

phenofile = rbind(file1, file2) %>%
  dplyr::mutate(names_in_RGset = paste0(.$Sentrix_ID, "_", .$Sentrix_Position),
                Basename = paste0(.$Sentrix_ID, "_", .$Sentrix_Position))


targets <- phenofile 

dat_daniela <- foreign::read.spss(file = "../spss epigenetica DS310121.sav", to.data.frame = T)

###information for the files
targets <- phenofile %>% 
  mutate(Sample_Group = ifelse(is.na(.$Sample_Group), .$Sample_Plate, .$Sample_Group),
         study = ifelse(grepl("MASHpool", .$Sample_Name), "MASHpool", ifelse(grepl("MASH", .$Sample_Name), "MASH", "BARIA")),
         sample_type = ifelse(grepl("subc", .$Sample_Name), "subc", 
                              ifelse(grepl("vis", .$Sample_Name), "vis",
                                     ifelse(grepl("liver", .$Sample_Name), "liver", 
                                            ifelse(Batch == "File2", "liver", "PBMC?")))),
         Subject_ID = ifelse(grepl("BARIA", .$Sample_Name), gsub("_.*","", gsub("BARIA_", "", gsub("replicate_","", .$Sample_Name))),
                             sub("_[^_]+$", "", .$Sample_Name)),
         her_OK = ifelse(grepl("her OK", .$Sample_Name), T, F)) %>%
  mutate(OK = ifelse(.$her_OK == T, F, ifelse(grepl("_OK", .$Sample_Name), T, F)),
         sex = dat_daniela$sex[match(.$Subject_ID, dat_daniela$ID)],
         age = dat_daniela$Age[match(.$Subject_ID, dat_daniela$ID)]) %>% 
  mutate(sex = ifelse(sex == "2" , "F", "M")) %>% 
  .[order(as.numeric(.$Subject_ID)),] %>% #gives the warning
  dplyr::select(-c(Pool_ID, Sample_Plate)) 

#Only use the BARIA cohort with liver DNA methylation
targets_baria <- targets[targets$study=="BARIA",] %>% 
  filter(sample_type == "liver") %>% 
  filter(!Sample_Name %in% c("BARIA_10_her OK", "BARIA_221_her OK", "BARIA_335_her OK")) %>% #Not in list
  filter(!Sample_Name == "BARIA_346_OK") #Sample swap (metadata listed as female with female medication; epigenetics data is of a male).

### load raw data from idats
RGset=minfi::read.metharray.exp(targets = targets_baria, idatspath, recursive = T)

pd <- pData(RGset) 
rownames(RGset@colData) <- pd$Sample_Name[match(rownames(RGset@colData), pd$names_in_RGset)]
x=getAnnotation(RGset)
length(rownames(x))

# Preprocessing
mset <- minfi::preprocessNoob(RGset)
gmset <- minfi::mapToGenome(mset)

#### Quality control of the data:
#Sex check
minfiQCqc <- minfi::minfiQC(gmset)
sex_dat <- minfiQCqc$qc %>% as.data.frame()

#age
beta_values <- data.frame(ID = colnames(gmset), t(minfi::getBeta(gmset)))

#### age prediction (error margin of 3-5 years); works for blood; might not work for liver tissue
bval_mat <- beta_values %>% filter(rownames(.) %in% targets_baria$Sample_Name) %>% dplyr::select(-ID) %>% as.matrix()
age_predicted <- agep(t(bval_mat), coeff = NULL, method = 'horvath')

t1 <- targets_baria %>%
  mutate(age_predicted = age_predicted$horvath.age[match(.$Sample_Name, rownames(age_predicted))],
         Predicted_sex = sex_dat$predictedSex[match(.$Sample_Name, rownames(sex_dat))],
         xmed = sex_dat$xMed[match(.$Sample_Name, rownames(sex_dat))],
         ymed = sex_dat$yMed[match(.$Sample_Name, rownames(sex_dat))]) %>% 
  mutate(label = ifelse(.$xmed<=6, .$Sample_Name, ""))
t1 %>% dplyr::select(Sample_Name, sex, Predicted_sex, age, age_predicted, sample_type) %>% View()

rmse = sqrt(mean((t1$age - t1$age_predicted)^2, na.rm=T))

path_qc_loc <- paste0(workdirectory, "../../Figures/QC")
dir.create(path_qc_loc)

## Age plot
age_plot <-
  ggscatter(data = t1,
            x = "age", y = "age_predicted",
            xlab = "Age",
            ylab = "Predicted age (Horvath)",
            size = 1) +
  stat_cor(method = "pearson",
           cor.coef.name = "rho",
           label.x = 40, label.y = 65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  annotate(geom = 'text',  label = paste("RMSE = ", round(rmse,2)),
           x= 40, y=67,
           hjust=0,
           fontface = "plain", color = "black")
ggsave(filename = "Age_prediction.png", plot = age_plot, path = path_qc_loc, device = "png")


## BEFROE REMOVAL OF STUFF (GAPS/SNPS/ALOSOMAL CpGs)
fin.pca <- prcomp(as.matrix(bval_mat), center = T, scale. = T)
s <- summary(fin.pca)

pc_percent <- s$importance[2,1:ncol(s$importance)]
pca_stuff <- as.data.frame(fin.pca$x) %>% 
  dplyr::select(1,2) %>% 
  setNames(c("X","Y")) %>% 
  dplyr::mutate(Sample = rownames(.)) %>% 
  dplyr::mutate(Batch = as.factor(phenofile$Batch[match(.$Sample, phenofile$Sample_Name)]),
                data = as.factor(targets$study[match(.$Sample, targets$Sample_Name)]),
                sample_type = as.factor(targets$sample_type[match(.$Sample, targets$Sample_Name)]))#,

#Check for batch effect (after only preprocessNoob) -> 1st component separates the batches.
p1_batch <- ggplot(pca_stuff, aes(x=X, y=Y))+
  geom_point(aes(color=Batch), size = 1)+
  xlab(paste0('PC1 [', round(pc_percent[1]*100,2),'%]'))+
  ylab(paste0('PC2 [', round(pc_percent[2]*100,2),'%]'))
p1_batch
ggsave(filename = "raw_pca_beta_values_after_preprocessNoob_batch", plot = p1_batch, device = "png", path = path_qc_loc)

#Sanity check for sex according to the methylation data and metadata
ps2a <-ggplot(t1, aes(x = xmed, y = ymed)) +
  geom_point(aes(color = sex), size = 1)+
  ggrepel::geom_label_repel(data = t1[t1$ymed>11 & t1$sex=="F",],
                            aes(x=xmed,y=ymed-0.3,label = Sample_Name)) +
  xlim(11, 14) +
  ylim(8, 14) +
  labs(title = "Sex data metadata", color='Sex')

ps2b <- ggplot(t1, aes(x = xmed, y = ymed)) +
  geom_point(aes(color = Predicted_sex), size = 1)+
  xlim(11, 14) +
  ylim(8, 14) +
  labs(title = "Sex data Predicted", color='Sex') +
  theme(legend.position = "none")

#Combine the plots
pl2 <- ps2a + ps2b+ plot_layout(guides="collect") + plot_annotation(tag_levels = c("A","B","C"))
pl2

#Add batch information to a third plot
ps3 <- ggplot(t1, aes(x = xmed, y = ymed)) +
  geom_point(aes(color = Batch), size = 1) +
  xlim(11, 14) +
  ylim(8, 14) +
  scale_color_manual(values=c("#CC6666", "#9999CC")) +
  labs(title = "Sample_type", color='Sample_type')
ps3 

p3 <- pl2 + ps3
ggsave(filename = "initial_batch_sexdata_plot_clean.png",plot = p3, device = "png",width = 10,height = 8, path = path_qc_loc)

#SNP sanity check on repeated measures
ewastools_object <- ewastools::read_idats(paste0(idatspath, targets_baria$Basename))
ewastools_betas <- ewastools::dont_normalize(ewastools_object)

pw_sample_corsnps <- cor(ewastools_betas[ewastools_object$manifest[probe_type == 'rs', index], ], use="complete.obs")
colnames(pw_sample_corsnps) <- targets_baria$Sample_Name[match(targets_baria$names_in_RGset, colnames(pw_sample_corsnps))]
pl <- pheatmap::pheatmap(pw_sample_corsnps)

pdf(file = paste0(path_qc_loc, "/heatmap_SNPs_baria_full.pdf"), width = 20, height = 20)
pl
dev.off()

#all duplicated_samples (HER_OK samples)
dup <- targets_baria[targets_baria$Subject_ID %in% targets_baria[duplicated(targets_baria$Subject_ID),]$Subject_ID,]
dup_plot <- pw_sample_corsnps[rownames(pw_sample_corsnps) %in% dup$Basename,
                              colnames(pw_sample_corsnps) %in% dup$Sample_Name]

pdf(file =  paste0(path_qc_loc, "/heatmap_SNPs_baria_repeated_measures.pdf"))
pheatmap::pheatmap(dup_plot)
dev.off()

## some processing parameters
gap_threshold = 0.2
min_maf = 0.01
snp_removal=T
gap_removal=T
allosome_removal=T

# Remove SNP probes
if(snp_removal){
  snp_cpgs <- minfi::getSnpInfo(gmset)[which(minfi::getSnpInfo(gmset)$CpG_maf > min_maf | minfi::getSnpInfo(gmset)$SBE_maf > min_maf),]
  gmset <- gmset[!rownames(gmset) %in% rownames(snp_cpgs),]
}

# Remove gapped probes
if(gap_removal){
  gmset_gh <- minfi::gaphunter(gmset, threshold = gap_threshold)
  gmset <- gmset[!rownames(gmset) %in% rownames(gmset_gh$sampleresults),]
}

# Annotations
annotations <- minfi::getAnnotation(gmset)

# Remove allosomal probes
if(allosome_removal){
  non_allosomal_probes <- rownames(annotations)[!annotations$chr %in% c("chrX", "chrY")]
  gmset <- gmset[rownames(gmset) %in% non_allosomal_probes,]
  annotations <- annotations[rownames(annotations) %in% non_allosomal_probes,]
}

genenames <- unlist(lapply(lapply(strsplit(annotations$UCSC_RefGene_Name, ";"), unique),
                           function(genes){
                             paste(genes, collapse = ".")
                           })
)

featurenames <- unlist(lapply(lapply(strsplit(annotations$UCSC_RefGene_Group, ";"), unique),
                              function(gfeats){
                                paste(gfeats, collapse = ".")
                              })
)

coordinates <- paste0(annotations$chr, ".", annotations$pos)

# CpG annotations
cpg_annotation <- paste0(rownames(gmset), "_", coordinates, "_", genenames, "_", featurenames)
cpg_annotation <- gsub("(^.+?)_+$", "\\1", cpg_annotation)

# Combat input
beta_values <- data.frame(ID = colnames(gmset), t(minfi::getBeta(gmset)))
colnames(beta_values) <- c("ID", cpg_annotation)
bval_mat <- beta_values %>% dplyr::select(-ID) %>% as.matrix()

# Combat normalization (functional normalization did not do the trick)
batch = targets_baria %>% column_to_rownames("Sample_Name") %>% dplyr::select(Batch) %>% .$Batch
mod = model.matrix(~1, data=targets_baria)
combat_edata = sva::ComBat(dat=t(bval_mat), batch=batch, mod=mod, par.prior=TRUE, prior.plots=T)

# Save beta-values (combat normalized)
dat <- list(X = t(combat_edata))[[1]]
rio::export(dat, paste0(workdirectory, "../../Intermediate_files/df_beta_val_noob_combat1.csv"), row.names=T)

### Same procedure for M-values
# re-done combat -> if I convert using log2(dat / (1 - dat)), NAs are introduced; PCA is very similar.
m_values <- data.frame(ID = colnames(gmset), t(minfi::getM(gmset)))
colnames(m_values) <- c("ID", cpg_annotation)
mval_mat <- m_values %>% dplyr::select(-ID) %>% as.matrix()
combat_edata_m <- sva::ComBat(dat=t(mval_mat), batch=batch, mod=mod, par.prior=TRUE, prior.plots=T)

# Save M-values (combat normalized)
dat_m <- list(X = t(combat_edata_m))[[1]]
rio::export(dat_m, paste0(workdirectory, "../../Intermediate_files/df_m_val_noob_combat.csv"), row.names=T)


#Pca after combat (beta-values)
fin.pca <- prcomp(as.matrix(t(combat_edata)), center = T, scale. = T)
plot(fin.pca, type = "l")
s <- summary(fin.pca)

pc_percent <- s$importance[2,1:ncol(s$importance)]
pca_stuff <- as.data.frame(fin.pca$x) %>% 
  dplyr::select(1,2) %>% 
  setNames(c("X","Y")) %>% 
  dplyr::mutate(Sample = rownames(.)) %>% 
  dplyr::mutate(Batch = as.factor(phenofile$Batch[match(.$Sample, phenofile$Sample_Name)]),
                data = as.factor(targets$study[match(.$Sample, targets$Sample_Name)]),
                sample_type = as.factor(targets$sample_type[match(.$Sample, targets$Sample_Name)]))
p1_combat <- ggplot(pca_stuff, aes(x=X, y=Y))+
  geom_point(aes(color=Batch), size = 1)+
  xlab(paste0('PC1 [', round(pc_percent[1]*100,2),'%]'))+
  ylab(paste0('PC2 [', round(pc_percent[2]*100,2),'%]'))

p1_combat
ggsave(filename = "PCA_Baria_ComBat.png", plot = p1_combat, device = "png", path = path_qc_loc)

### pca after combat (M-values)
fin.pca <- prcomp(as.matrix(t(combat_edata_m)), center = T, scale. = T)
plot(fin.pca, type = "l")
s <- summary(fin.pca)

pc_percent <- s$importance[2,1:ncol(s$importance)]
pca_stuff <- as.data.frame(fin.pca$x) %>% 
  dplyr::select(1,2) %>% 
  setNames(c("X","Y")) %>% 
  dplyr::mutate(Sample = rownames(.)) %>% 
  dplyr::mutate(Batch = as.factor(phenofile$Batch[match(.$Sample, phenofile$Sample_Name)]),
                data = as.factor(targets$study[match(.$Sample, targets$Sample_Name)]),
                sample_type = as.factor(targets$sample_type[match(.$Sample, targets$Sample_Name)]))
p1_combat_m <- ggplot(pca_stuff, aes(x=X, y=Y))+
  geom_point(aes(color=Batch), size = 1)+
  xlab(paste0('PC1 [', round(pc_percent[1]*100,2),'%]'))+
  ylab(paste0('PC2 [', round(pc_percent[2]*100,2),'%]'))

p1_combat_m
ggsave(filename = "PCA_Baria_ComBat_m.png", plot = p1_combat_m, device = "png", path = path_qc_loc)