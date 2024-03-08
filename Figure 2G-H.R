###codes for Figure 2G-H.

#1 Figure 2G.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
library(Seurat)
HNSCC <- readRDS(file = "HNSCC_resolution0.3_annotation_sub.rds") #The detailed code for HNSCC_resolution0.3_annotation_sub.rds is at line 273 of the Figure 1C-F.R file.
table(HNSCC$Celltype)
table(HNSCC$Sample)

###Just tumor samples.
HNSCC_T <- subset(HNSCC, subset=Sample %in% c("Primary", "Metastatic")) 
HNSCC_T$Sample <- as.character(as.factor(HNSCC_T$Sample))
HNSCC_T$orig.ident <- as.character(as.factor(HNSCC_T$orig.ident))
table(HNSCC_T$Sample)
table(HNSCC_T$orig.ident)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/HNSCC")
gs = lapply(readLines("HNSCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
HNSCC_T <- AddModuleScore(
  object = HNSCC_T,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(HNSCC_T@meta.data)[9] <- 'apCAFs'
colnames(HNSCC_T@meta.data)[10] <- 'CD4T'
head(HNSCC_T@meta.data)

###Retain only the necessary data.
m <- HNSCC_T@meta.data
df_HNSCC_T <- m[, c(9, 10)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_HNSCC_T$apCAFs)
ad.test(df_HNSCC_T$CD4T)

###Visualization.
#scatter plot for Figure 2G.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
pdf(file = "HNSCC_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_HNSCC_T, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Figure 2H.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
library(Seurat)
OV <- readRDS(file = "OV_resolution0.2_annotation2.rds") #The detailed code for OV_resolution0.2_annotation2.rds is at line 261 of the Figure 1G-J.R file.
table(OV$Celltype)
table(OV$Sample)
table(OV$Treatment)

###Just OV tumors samples without drug treatment.
OV_naive <- subset(OV, subset=Treatment %in% c("naive")) 
OV_naive$Sample <- as.character(as.factor(OV_naive$Sample))
OV_naive$orig.ident <- as.character(as.factor(OV_naive$orig.ident))
OV_naive$Treatment <- as.character(as.factor(OV_naive$Treatment))
table(OV_naive$Sample)
table(OV_naive$orig.ident)
table(OV_naive$Treatment)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/OV")
gs = lapply(readLines("OV_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
OV_naive <- AddModuleScore(
  object = OV_naive,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
colnames(OV_naive@meta.data)[10] <- 'apCAFs'
colnames(OV_naive@meta.data)[11] <- 'CD4T'
head(OV_naive@meta.data)

###Retain only the necessary data.
m <- OV_naive@meta.data
df_OV_naive <- m[, c(10, 11)]

###Correlation analysis.
###Anderson-Darling test for normality.
library(nortest)
library(ggpubr)
ad.test(df_OV_naive$apCAFs)
ad.test(df_OV_naive$CD4T)

###Visualization.
#scatter plot for Figure 2H.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
pdf(file = "OV_SC_apCAFs_effector_CD4_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_OV_naive, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
