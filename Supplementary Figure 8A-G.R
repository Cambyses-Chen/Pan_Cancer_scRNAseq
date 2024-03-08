###codes for Supplementary Figure 8A-G.

#1 Supplementary Figure 8A.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
library(Seurat)
NPC <- readRDS(file = "NPC_resolution0.2_annotation_sub.rds") #The detailed code for NPC_resolution0.2_annotation_sub.rds is at line 616 of the Supplementary Figure 5A-B.R file.
table(NPC$Celltype)
table(NPC$Sample)

###Just tumor samples.
NPC_T <- subset(NPC, subset=Sample %in% c("Tumor")) 
NPC_T$Sample <- as.character(as.factor(NPC_T$Sample))
NPC_T$orig.ident <- as.character(as.factor(NPC_T$orig.ident))
table(NPC_T$Sample)
table(NPC_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/NPC")
gs = lapply(readLines("NPC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
NPC_T <- AddModuleScore(
  object = NPC_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(NPC_T@meta.data)[9] <- 'apCAFs'
colnames(NPC_T@meta.data)[10] <- 'CD4T'
head(NPC_T@meta.data)

###Retain only the necessary data.
m <- NPC_T@meta.data
df_NPC_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_NPC_T$apCAFs)
ad.test(df_NPC_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8A.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
pdf(file = "NPC_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_NPC_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 8B.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
library(Seurat)
BRCA <- readRDS(file = "BC_GSE176078_resolution0.1_annotation_sub.rds") #The detailed code for BC_GSE176078_resolution0.1_annotation_sub.rds is at line 308 of the Supplementary Figure 5G-H.R file.
table(BRCA$Celltype)
table(BRCA$Sample)

###Just tumor samples.
BRCA_T <- subset(BRCA, subset=Sample %in% c("ER", "HER2", "TNBC")) 
BRCA_T$Sample <- as.character(as.factor(BRCA_T$Sample))
BRCA_T$orig.ident <- as.character(as.factor(BRCA_T$orig.ident))
table(BRCA_T$Sample)
table(BRCA_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/BRCA")
gs = lapply(readLines("BRCA_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
BRCA_T <- AddModuleScore(
  object = BRCA_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(BRCA_T@meta.data)[9] <- 'apCAFs'
colnames(BRCA_T@meta.data)[10] <- 'CD4T'
head(BRCA_T@meta.data)

###Retain only the necessary data.
m <- BRCA_T@meta.data
df_BRCA_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_BRCA_T$apCAFs)
ad.test(df_BRCA_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8B.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
pdf(file = "BRCA_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_BRCA_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#3 Supplementary Figure 8C.
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
library(Seurat)
AM <- readRDS(file = "AM_resolution0.1_annotation_sub.rds") #The detailed code for AM_resolution0.1_annotation_sub.rds is at line 436 of the Supplementary Figure 5I-J.R file.
table(AM$Celltype)
table(AM$Sample)

###Just tumor samples.
AM_T <- subset(AM, subset=Sample %in% c("primary", "metastatic")) 
AM_T$Sample <- as.character(as.factor(AM_T$Sample))
AM_T$orig.ident <- as.character(as.factor(AM_T$orig.ident))
table(AM_T$Sample)
table(AM_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/AM")
gs = lapply(readLines("AM_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
AM_T <- AddModuleScore(
  object = AM_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(AM_T@meta.data)[9] <- 'apCAFs'
colnames(AM_T@meta.data)[10] <- 'CD4T'
head(AM_T@meta.data)

###Retain only the necessary data.
m <- AM_T@meta.data
df_AM_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_AM_T$apCAFs)
ad.test(df_AM_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8C.
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
pdf(file = "AM_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_AM_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#4 Supplementary Figure 8D.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
library(Seurat)
CM <- readRDS(file = "CM_resolution0.2_annotation2.rds") #The detailed code for CM_resolution0.2_annotation2.rds is at line 246 of the Supplementary Figure 5K-L.R file.
table(CM$Celltype)
table(CM$Sample)

###Just tumor samples.
CM_T <- subset(CM, subset=Sample %in% c("primary", "metastatic")) 
CM_T$Sample <- as.character(as.factor(CM_T$Sample))
CM_T$orig.ident <- as.character(as.factor(CM_T$orig.ident))
table(CM_T$Sample)
table(CM_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/CM")
gs = lapply(readLines("CM_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
CM_T <- AddModuleScore(
  object = CM_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(CM_T@meta.data)[9] <- 'apCAFs'
colnames(CM_T@meta.data)[10] <- 'CD4T'
head(CM_T@meta.data)

###Retain only the necessary data.
m <- CM_T@meta.data
df_CM_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_CM_T$apCAFs)
ad.test(df_CM_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8D.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
pdf(file = "CM_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_CM_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#5 Supplementary Figure 8E.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
library(Seurat)
RCC <- readRDS(file = "RCC_resolution0.1_annotation2.rds") #The detailed code for RCC_resolution0.1_annotation2.rds is at line 288 of the Supplementary Figure 5M-N.R file.
table(RCC$Celltype)
table(RCC$Sample)

###Just tumor samples.
RCC_T <- subset(RCC, subset=Sample %in% c("Primary tumor")) 
RCC_T$Sample <- as.character(as.factor(RCC_T$Sample))
RCC_T$orig.ident <- as.character(as.factor(RCC_T$orig.ident))
table(RCC_T$Sample)
table(RCC_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/RCC")
gs = lapply(readLines("RCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
RCC_T <- AddModuleScore(
  object = RCC_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(RCC_T@meta.data)[9] <- 'apCAFs'
colnames(RCC_T@meta.data)[10] <- 'CD4T'
head(RCC_T@meta.data)

###Retain only the necessary data.
m <- RCC_T@meta.data
df_RCC_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_RCC_T$apCAFs)
ad.test(df_RCC_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8E.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
pdf(file = "RCC_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_RCC_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#6 Supplementary Figure 8F.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
library(Seurat)
CC <- readRDS(file = "Cervical_resolution0.3_annotation_2.rds") #The detailed code for Cervical_resolution0.3_annotation_2.rds is at line 268 of the Supplementary Figure 5C-D.R file.
table(CC$Celltype)
table(CC$Sample)

###Just tumor samples.
CC_T <- subset(CC, subset=Sample %in% c("Tumor")) 
CC_T$Sample <- as.character(as.factor(CC_T$Sample))
CC_T$orig.ident <- as.character(as.factor(CC_T$orig.ident))
table(CC_T$Sample)
table(CC_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/CC")
gs = lapply(readLines("CC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
CC_T <- AddModuleScore(
  object = CC_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(CC_T@meta.data)[9] <- 'apCAFs'
colnames(CC_T@meta.data)[10] <- 'CD4T'
head(CC_T@meta.data)

###Retain only the necessary data.
m <- CC_T@meta.data
df_CC_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_CC_T$apCAFs)
ad.test(df_CC_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8F.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
pdf(file = "CC_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_CC_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#7 Supplementary Figure 8G.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
library(Seurat)
CRC <- readRDS(file = "CRC_HRA000979_resolution0.1_annotation_sub.rds") #The detailed code for CRC_HRA000979_resolution0.1_annotation_sub.rds is at line 269 of the Supplementary Figure 5E-F.R file.
table(CRC$Celltype)
table(CRC$Sample)

###Just tumor samples.
CRC_T <- subset(CRC, subset=Sample %in% c("Tumor")) 
CRC_T$Sample <- as.character(as.factor(CRC_T$Sample))
CRC_T$orig.ident <- as.character(as.factor(CRC_T$orig.ident))
table(CRC_T$Sample)
table(CRC_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/CRC")
gs = lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
CRC_T <- AddModuleScore(
  object = CRC_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(CRC_T@meta.data)[9] <- 'apCAFs'
colnames(CRC_T@meta.data)[10] <- 'CD4T'
head(CRC_T@meta.data)

###Retain only the necessary data.
m <- CRC_T@meta.data
df_CRC_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_CRC_T$apCAFs)
ad.test(df_CRC_T$CD4T)

###Visualization.
#scatter plot for Supplementary Figure 8G.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
pdf(file = "CRC_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_CRC_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
