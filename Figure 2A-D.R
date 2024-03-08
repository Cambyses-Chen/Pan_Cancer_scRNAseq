###codes for Figure 2A-D.

#1 Figure 2A, C.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
library(Seurat)
library(tidyverse)
HNSCC = readRDS(file = "HNSCC_Tumor_merge.rds") ##HNSCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of HNSCC tumors. The specific code is in the file named HNSCC_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(HNSCC) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(HNSCC,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

###The gene of p_val_adj < 0.01 was retained.
ref_marker=ref_marker%>%filter(p_val_adj < 0.01)

###The gene of d > 0.2 was retained.
ref_marker$d=ref_marker$pct.1 - ref_marker$pct.2
ref_marker=ref_marker%>%filter(d > 0.2)

###The genes are sequenced in descending order by avg_log2FC.
ref_marker=ref_marker%>%arrange(cluster,desc(avg_log2FC))
ref_marker=as.data.frame(ref_marker)

###Extraction of signature genes.
used_gene=sort(unique(ref_marker$gene))

###Extraction of signature matrix.
sig_matrix=AverageExpression(HNSCC, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "HNSCC_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-HNSC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-HNSC")
library(TCGAbiolinks)
load("TCGA-HNSC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "tpm_unstrand")
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
saveRDS(t10, file = "TCGA-HNSC_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-HNSC_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-HNSC_normal_RNAseq_TPM.rds")

###The genes in the TCGA-HNSC tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(HNSCC@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-HNSC Mixture.
saveRDS(t13, file = "TCGA-HNSC_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-HNSC Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-HNSC_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-HNSC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-HNSC_Results.csv")

###Loading TCGA-HNSC clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-HNSC")
clin2 <- readRDS(file = "TCGA-HNSC_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

###Data preprocessing.
t14 <- t(t13)
t14 <- as.data.frame(t14)
t14$id <- rownames(t14)

#Remove duplicate id.
t15 <- t14[!duplicated(t14$id),]

#The samples in the clinical data must match those in the expression matrix.
library(stringr)
t16 <- clin2[clin2[,1] %in% str_sub(rownames(t15),1,12),]

#Remove duplicate bcr_patient_barcode.
t17 <- t16[!duplicated(t16$bcr_patient_barcode),]

#Retain only the necessary clinical information.
t17 <- t17[, c(1, 8, 11, 12)]

#Extract clinical information for alive and dead patients separately.
t18 <- t17[grep("Alive", t17[,2]),]
t19 <- t17[grep("Dead", t17[,2]),]
t18$OS.time <- t18$days_to_last_followup
t19$OS.time <- t19$days_to_death
t18$status <- 0
t19$status <- 1
t20 <- rbind(t18, t19)
t20$barcode <- t20$bcr_patient_barcode
t20 <- t20[, c(7, 6, 5)]
t21 <- t20[!is.na(t20$OS.time),]

#Merge sample clinical information with sample expression matrix.
t15$barcode <- str_sub(rownames(t15),1,12)
t22 <- t15[!duplicated(t15$barcode),]
t23 <- t22[t22$barcode %in% t21$barcode,]
t24 <- merge(t21, t23, by.x = "barcode", by.y = "barcode")

#Convert barcode to id.
rownames(t24) <- t24$id
t24$barcode <- t24$id
saveRDS(t24, file = "TCGA-HNSC_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-HNSC_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:5), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:47),]
s5$group <- "low"
s6<- s4[c(48:514),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Figure 2A.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 467','Low  n = 47'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-HNSC_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
print(survp, newpage = FALSE)
dev.off()

##5 Correlation analysis.
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor.cells)

###Visualization.
#scatter plot for Figure 2C.
pdf(file = "TCGA-HNSC_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Figure 2B, D.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
library(Seurat)
library(tidyverse)
OV = readRDS(file = "OV_naive_merge.rds") ##OV_naive_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of OV tumors without drug treatment. The specific code is in the file named OV_naive_merge.R.

###Change the default cell identity to "celltype".
Idents(OV) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(OV,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

###The gene of p_val_adj < 0.01 was retained.
ref_marker=ref_marker%>%filter(p_val_adj < 0.01)

###The gene of d > 0.2 was retained.
ref_marker$d=ref_marker$pct.1 - ref_marker$pct.2
ref_marker=ref_marker%>%filter(d > 0.2)

###The genes are sequenced in descending order by avg_log2FC.
ref_marker=ref_marker%>%arrange(cluster,desc(avg_log2FC))
ref_marker=as.data.frame(ref_marker)

###Extraction of signature genes.
used_gene=sort(unique(ref_marker$gene))

###Extraction of signature matrix.
sig_matrix=AverageExpression(OV, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "OV_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-OV RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-OV")
load("TCGA-OV_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "tpm_unstrand")
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
saveRDS(t10, file = "TCGA-OV_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-OV_tumor_RNAseq_TPM.rds")

###The genes in the TCGA-OV tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(OV@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-OV Mixture.
saveRDS(t13, file = "TCGA-OV_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-OV Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-OV_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-OV")
ciber <- read.csv(file = "CIBERSORTx_TCGA-OV_Results.csv")

###Loading TCGA-OV clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-OV")
clin2 <- readRDS(file = "TCGA-OV_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

###Data preprocessing.
t14 <- t(t13)
t14 <- as.data.frame(t14)
t14$id <- rownames(t14)

#Remove duplicate id.
t15 <- t14[!duplicated(t14$id),]

#The samples in the clinical data must match those in the expression matrix.
library(stringr)
t16 <- clin2[clin2[,1] %in% str_sub(rownames(t15),1,12),]

#Remove duplicate bcr_patient_barcode.
t17 <- t16[!duplicated(t16$bcr_patient_barcode),]

#Retain only the necessary clinical information.
t17 <- t17[, c(1, 6, 8, 9)]

#Extract clinical information for alive and dead patients separately.
t18 <- t17[grep("Alive", t17[,2]),]
t19 <- t17[grep("Dead", t17[,2]),]
t18$OS.time <- t18$days_to_last_followup
t19$OS.time <- t19$days_to_death
t18$status <- 0
t19$status <- 1
t20 <- rbind(t18, t19)
t20$barcode <- t20$bcr_patient_barcode
t20 <- t20[, c(7, 6, 5)]
t21 <- t20[!is.na(t20$OS.time),]

#Merge sample clinical information with sample expression matrix.
t15$barcode <- str_sub(rownames(t15),1,12)
t22 <- t15[!duplicated(t15$barcode),]
t23 <- t22[t22$barcode %in% t21$barcode,]
t24 <- merge(t21, t23, by.x = "barcode", by.y = "barcode")

#Convert barcode to id.
rownames(t24) <- t24$id
t24$barcode <- t24$id
saveRDS(t24, file = "TCGA-OV_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-OV_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]

###Arrange in increasing order according to apCAFs scores.
s4 <- s2[order(s2$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:152),]
s5$group <- "low"
s6<- s4[c(153:419),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Figure 2B.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 267','Low  n = 152'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-OV_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
print(survp, newpage = FALSE)
dev.off()

##5 Correlation analysis.
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor.cells)

###Visualization.
#scatter plot for Figure 2D.
pdf(file = "TCGA-OV_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
