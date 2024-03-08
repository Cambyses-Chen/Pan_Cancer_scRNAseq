###codes for Figure 2E-F.
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)

#1 Figure 2E.
###Loading and processing TCGA-HNSC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-HNSC")
library(TCGAbiolinks)
load("TCGA-HNSC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 568)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-HNSC_RNAseq_count.rds")

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
saveRDS(t11, file = "TCGA-HNSC_tumor_RNAseq_count.rds")
saveRDS(t12, file = "TCGA-HNSC_normal_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature. 
setwd("E:/Raw Data/Signatures/HNSCC")
gs = lapply(readLines("HNSCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("HNSCC_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "HNSCC_TCGA-HNSC_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Figure 2E.
pdf(file = "HNSCC_TCGA-HNSC_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Figure 2F.
###Loading and processing TCGA-OV RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-OV")
load("TCGA-OV_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 431)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-OV_RNAseq_count.rds")

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
t11 <- t10[,tumor_id]
t11 <- as.matrix(t11)

###Tumor sample data were stored.
saveRDS(t11, file = "TCGA-OV_tumor_RNAseq_count.rds")

data <- t11

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/OV")
gs = lapply(readLines("OV_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Run GSVA.
es.dif <- gsva(data, gs, mx.diff=TRUE, kcdf="Poisson",method = "gsva", verbose=FALSE, parallel.sz=10)

###Offer the rowname.
rownames(es.dif) = unlist(lapply(readLines("OV_apCAFs_CD4.txt"), function(x){
  strsplit(x, "\t")[[1]][1]
}))
saveRDS(es.dif, file = "OV_TCGA-OV_apCAFs_CD4_effector_T_cells_es.dif.rds")

###Correlation analysis.
data2 <- t(es.dif)
data2 <- as.data.frame(data2)

###Shapiro-Wilk test for normality.
shapiro.test(data2$apCAFs)
shapiro.test(data2$CD4_T)

###Visualization.
#scatter plot for Figure 2F.
pdf(file = "OV_TCGA-OV_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(data2, x = "CD4_T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
