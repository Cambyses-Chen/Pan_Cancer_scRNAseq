###codes for Supplementary Figure 7O-U.
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)

#1 Supplementary Figure 7O.
###Load counts matrix.
setwd("E:/Raw Data/TCGA_data/GSE102349_NPC_RNAseq")
data <- readRDS(file = "GSE102349_NPC_RNAseq.rds") #Downloaded and processed from the GEO database.
data <- as.matrix(data)

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/NPC")
gs = lapply(readLines("NPC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Gaussian",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("NPC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "NPC_GSE102349_NPC_RNAseq_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7O.
pdf(file = "GSE102349_NPC_RNAseq_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 7P.
###Loading and processing TCGA-CESC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-CESC")
load("TCGA-CESC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[, c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 311)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-CESC_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,14) == "0","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-CESC_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-CESC_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/CC")
gs = lapply(readLines("CC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("CC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "CC_TCGA-CESC_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7P.
pdf(file = "CC_TCGA-CESC_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#3 Supplementary Figure 7Q.
###Loading and processing TCGA-COAD RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-COAD")
load("TCGA-COAD_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[, c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 526)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-COAD_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,14) == "0","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-COAD_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-COAD_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/CRC")
gs = lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "CRC_TCGA-COAD_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data <- t(es.dif)
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7Q.
pdf(file = "CRC_TCGA-COAD_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#4 Supplementary Figure 7R.
setwd("E:/Raw Data/TCGA_data/TCGA-READ")
load("TCGA-READ_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[, c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 179)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-READ_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,14) == "0","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),                    
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-READ_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-READ_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/CRC")
gs = lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "CRC_TCGA-READ_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7R.
pdf(file = "CRC_TCGA-READ_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#5 Supplementary Figure 7S.
###Loading and processing TCGA-BRCA RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-BRCA")
load("TCGA-BRCA_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[,c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 1233)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-BRCA_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,15) == "01","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-BRCA_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-BRCA_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/BRCA")
gs = lapply(readLines("BRCA_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("BRCA_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "BRCA_TCGA-BRCA_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7S.
pdf(file = "BRCA_TCGA-BRCA_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#6 Supplementary Figure 7T.
###Loading and processing TCGA-SKCM RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-SKCM")
load("TCGA-SKCM_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[, c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 475)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-SKCM_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,14) == "0","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-SKCM_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-SKCM_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/CM")
gs = lapply(readLines("CM_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("CM_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "CM_TCGA-SKCM_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7T.
pdf(file = "CM_TCGA-SKCM_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#7 Supplementary Figure 7U.
###Loading and processing TCGA-KIRC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-KIRC")
load("TCGA-KIRC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "unstranded") #count data.
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[, c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 616)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-KIRC_RNAseq_count.rds")

###Distinguish between tumor and paracancer samples.
t10 <- as.data.frame(t10)
library(stringr)
group_list<-ifelse(str_sub(colnames(t10),14,14) == "0","tumor","normal")
group_list
coldata<-data.frame(row.names = colnames(t10),
                    condition=group_list)
coldata
coldata$id <- rownames(coldata)
tumor <- coldata[grep("tumor", coldata[,1]),]
normal <- coldata[grep("normal", coldata[,1]),]
tumor_id <- tumor$id
normal_id <- normal$id
t11 <- t10[,tumor_id]
t12 <- t10[,normal_id]
t11 <- as.matrix(t11)
t12 <- as.matrix(t12)

###Tumor and paracancer sample data were stored separately.
saveRDS(t11, file = "TCGA-KIRC_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-KIRC_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/RCC")
gs = lapply(readLines("RCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("RCC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "RCC_TCGA-KIRC_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Supplementary Figure 7U.
pdf(file = "RCC_TCGA-KIRC_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", #spearman#pearson
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
