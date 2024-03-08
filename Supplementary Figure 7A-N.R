###codes for Supplementary Figure 7A-N.

#1 Supplementary Figure 7A, H.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
library(Seurat)
library(tidyverse)
NPC = readRDS(file = "NPC_Tumor_merge.rds") #NPC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of NPC tumors. The specific code is in the file named NPC_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(NPC) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(NPC,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

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
sig_matrix=AverageExpression(NPC, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "NPC_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing GSE102349_NPC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/GSE102349_NPC_RNAseq")
t11 <- readRDS(file = "GSE102349_NPC_RNAseq.rds") #Downloaded and processed from the GEO database.
t11 <- as.data.frame(t11)
t11$gene_name <- rownames(t11)
t12 <- t11[!duplicated(t11$gene_name),]

###The genes in the GSE102349_NPC tumor sample matrix aligned with the genes in the Seurat object matrix.
Gene <- rownames(NPC@assays$RNA@meta.features)
t13 <- t12[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]
t13 <- t13[!duplicated(t13$gene_name),]
t13 <- t13[, -114]

###Save the GSE102349_NPC Mixture.
saveRDS(t13, file = "GSE102349_NPC_RNAseq_tumor_RNAseq_filterd.rds")

###Save the GSE102349_NPC Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "GSE102349-NPC_RNAseq_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/GSE102349_NPC_RNAseq")
ciber <- read.csv(file = "CIBERSORTx_GSE102349_NPC_RNAseq_Results.csv")

###Loading GSE102349_NPC clinical data.
setwd("E:/Raw Data/TCGA_data/GSE102349_NPC_RNAseq")
clin2 = read.table(file = "GSE102349.txt", header=T, sep="\t", check.names=F) #Downloaded and processed from the GEO database.

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% clin2$id,]
t26 <- merge(clin2, t25, by.x = "id", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "GSE102349_NPC_RNAseq_clinical_processed_cibersort_score.rds")

###Processing integrated data.
#Retain only the necessary clinical information. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

#Change names.
colnames(t27) <- c("id", "RFS.time", "status", "apCAFs")

#Convert months into days.
t27$RFS.time <- t27$RFS.time*30
s1 <- t27
s2 <- s1[order(s1$RFS.time),]

###Arrange in increasing order according to apCAFs scores.
s4 <- s2[order(s2$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:11),]
s5$group <- "low"
s6<- s4[c(12:88),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(RFS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7A.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 77','Low  n = 11'), legend = c(0.75, 0.2), risk.table = F, xlab='days', ylab='Progression Free Survival', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_GSE102349_NPC_RNAseq_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
#scatter plot for Supplementary Figure 7H.
pdf(file = "GSE102349_NPC_RNAseq_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 7B, I.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
library(Seurat)
library(tidyverse)
CC = readRDS(file = "Cervical_Tumor_merge.rds") #Cervical_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of Cervical cancers. The specific code is in the file named Cervical_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(CC) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(CC,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

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
sig_matrix=AverageExpression(CC, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "CC_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-CESC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-CESC")
load("TCGA-CESC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 311)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-CESC_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-CESC_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-CESC_normal_RNAseq_TPM.rds")

###The genes in the TCGA-CESC tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(CC@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-CESC Mixture.
saveRDS(t13, file = "TCGA-CESC_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-CESC Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-CESC_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-CESC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-CESC_Results.csv")

###Loading TCGA-CESC clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-CESC")
clin2 <- readRDS(file = "TCGA-CESC_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 21, 22, 23)]

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
saveRDS(t24, file = "TCGA-CESC_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-CESC_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:18), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:75),]
s5$group <- "low"
s6<- s4[c(76:286),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7B.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 211','Low  n = 75'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-CESC_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
#scatter plot for Supplementary Figure 7I.
pdf(file = "TCGA-CESC_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#3 Supplementary Figure 7C, D, J, K.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
library(Seurat)
library(tidyverse)
CRC = readRDS(file = "CRC_Tumor_merge.rds") #CRC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of CRC tumors. The specific code is in the file named CRC_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(CRC) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(CRC,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

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
sig_matrix=AverageExpression(CRC, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "CRC_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
#TCGA-COAD
###Loading and processing TCGA-COAD RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-COAD")
load("TCGA-COAD_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 524)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-COAD_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-COAD_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-COAD_normal_RNAseq_TPM.rds")

###The genes in the TCGA-COAD tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(CRC@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-COAD Mixture.
saveRDS(t13, file = "TCGA-COAD_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-COAD Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-COAD_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-COAD")
ciber <- read.csv(file = "CIBERSORTx_TCGA-COAD_Results.csv")

###Loading TCGA-COAD clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-COAD")
clin2 <- readRDS(file = "TCGA-COAD_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 7, 10, 11)]

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
saveRDS(t24, file = "TCGA-COAD_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-COAD_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:111), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:171),]
s5$group <- "low"
s6<- s4[c(172:343),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7C.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 172','Low  n = 171'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-COAD_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
#scatter plot for Supplementary Figure 7J.
pdf(file = "TCGA-COAD_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

##6 Making Mixture.
#TCGA-READ
###Loading and processing TCGA-READ RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-READ")
load("TCGA-READ_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 179)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-READ_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-READ_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-READ_normal_RNAseq_TPM.rds")

###The genes in the TCGA-READ tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(CRC@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-READ Mixture.
saveRDS(t13, file = "TCGA-READ_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-READ Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-READ_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##7 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##8 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-READ")
ciber <- read.csv(file = "CIBERSORTx_TCGA-READ_Results.csv")

###Loading TCGA-READ clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-READ")
clin2 <- readRDS(file = "TCGA-READ_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 7, 10, 11)]

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
saveRDS(t24, file = "TCGA-READ_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-READ_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:42), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:52),]
s5$group <- "low"
s6<- s4[c(53:123),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7D.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 71','Low  n = 52'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-READ_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
print(survp, newpage = FALSE)
dev.off()

##9 Correlation analysis.
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor.cells)

###Visualization.
#scatter plot for Supplementary Figure 7K.
pdf(file = "TCGA-READ_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#4 Supplementary Figure 7E, L.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
library(Seurat)
library(tidyverse)
BRCA = readRDS(file = "BC_Tumor_merge.rds") #BC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of BRCA tumors. The specific code is in the file named BC_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(BRCA) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(BRCA,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

###The gene of p_val_adj < 0.01 was retained.
ref_marker=ref_marker%>%filter(p_val_adj < 0.01)

###The gene of d > 0.2 was retained
ref_marker$d=ref_marker$pct.1 - ref_marker$pct.2
ref_marker=ref_marker%>%filter(d > 0.2)

###The genes are sequenced in descending order by avg_log2FC.
ref_marker=ref_marker%>%arrange(cluster,desc(avg_log2FC))
ref_marker=as.data.frame(ref_marker)

###Extraction of signature genes.
used_gene=sort(unique(ref_marker$gene))

###Extraction of signature matrix.
sig_matrix=AverageExpression(BRCA, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "BRCA_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-BRCA RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-BRCA")
load("TCGA-BRCA_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
library(TCGAbiolinks)
dataPrep <- TCGAanalyze_Preprocessing(object = data,
                                      cor.cut = 0.6, datatype = "tpm_unstrand")
t <- as.data.frame(dataPrep)
t$gene_id <- rownames(t)

###Replace the Ensemble ID with a SYMBOL ID.
t6 <- as.data.frame(data@rowRanges@elementMetadata@listData) 
t8 <- t6[,c(5, 7)]
t9 <- merge(t, t8, by.x = "gene_id", by.y = "gene_id")

###Remove duplicate SYMBOL ID.
t10 <- t9[!duplicated(t9$gene_name),]
rownames(t10) <- t10$gene_name
t10 <- t10[, -c(1, 1232)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-BRCA_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-BRCA_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-BRCA_normal_RNAseq_TPM.rds")

###The genes in the TCGA-BRCA tumor sample matrix aligned with the genes in the Seurat object matrix.
Gene <- rownames(BRCA@assays$RNA@meta.features)
t11 <- as.data.frame(t11)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-BRCA Mixture.
saveRDS(t13, file = "TCGA-BRCA_tumor_RNAseq_TPM_filterd.rds")

###Due to the large number of Mixture samples, they must be split into two separate parts.
p1 <- t13[, c(1:600)] 
p2 <- t13[, c(601:1117)]

###Store the Mixture of the two parts in txt format separately and upload them to the CIBERSORTx website for analysis.
write.table(p1,file = "TCGA-BRCA_p1_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)
write.table(p2,file = "TCGA-BRCA_p2_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)
write.table(t13,file = "TCGA-BRCA_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###Merge two parts of the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-BRCA")
ciber_p1 <- read.csv(file = "CIBERSORTx_TCGA-BRCA-p1_Results.csv")
ciber_p2 <- read.csv(file = "CIBERSORTx_TCGA-BRCA-p2_Results.csv")
ciber <- rbind(ciber_p1, ciber_p2)
saveRDS(ciber, file = "CIBERSORTx_TCGA-BRCA_Results.rds")

###Loading TCGA-BRCA clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-BRCA")
clin2 <- readRDS(file = "TCGA-BRCA_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 7, 10, 11)]

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
saveRDS(t24, file = "TCGA-BRCA_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-BRCA_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:61), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:323),]
s5$group <- "low"
s6<- s4[c(324:1033),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7E.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 710','Low  n = 323'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-BRCA_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
#scatter plot for Supplementary Figure 7L.
pdf(file = "TCGA-BRCA_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#5 Supplementary Figure 7F, M.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
library(Seurat)
library(tidyverse)
CM = readRDS(file = "CM_Tumor_merge.rds") #CM_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of CM tumors. The specific code is in the file named CM_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(CM) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(CM,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

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
sig_matrix=AverageExpression(CM, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "CM_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-SKCM RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-SKCM")
load("TCGA-SKCM_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 475)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-SKCM_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-SKCM_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-SKCM_normal_RNAseq_TPM.rds")

###The genes in the TCGA-SKCM tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(CM@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-SKCM Mixture.
saveRDS(t13, file = "TCGA-SKCM_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-SKCM Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-SKCM_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-SKCM")
ciber <- read.csv(file = "CIBERSORTx_TCGA-SKCM_Results.csv")

###Loading TCGA-SKCM clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-SKCM")
clin2 <- readRDS(file = "TCGA-SKCM_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 21, 22, 23)]

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
saveRDS(t24, file = "TCGA-SKCM_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-SKCM_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:24), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:20),]
s5$group <- "low"
s6<- s4[c(21:435),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7F.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 415','Low  n = 20'), legend = c(0.75, 0.85), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-SKCM_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
#scatter plot for Supplementary Figure 7M.
pdf(file = "TCGA-SKCM_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#6 Supplementary Figure 7G, N.
##1 Making signature matrix.
###Load Seurat object.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
library(Seurat)
library(tidyverse)
RCC = readRDS(file = "RCC_Tumor_merge.rds") #RCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of RCC tumors. The specific code is in the file named RCC_Tumor_merge.R.

###Change the default cell identity to "celltype".
Idents(RCC) = "celltype"

###Gene expression markers for all celltype.
ref_marker=FindAllMarkers(RCC,logfc.threshold = 0.8,min.pct = 0.1,only.pos = T,test.use = "wilcox")

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
sig_matrix=AverageExpression(RCC, slot = 'data' ,verbose = FALSE,features = used_gene)

###The multiplication by 100 here is necessary because Seurat's scale.factor is 10000, while TPM/CPM is defined as per one million. To maintain consistency, we multiply by 100 again.
sig_matrix=sig_matrix$RNA * 100
sig_matrix=as.data.frame(sig_matrix)

###Save the signature matrix in txt format and upload it to the CIBERSORTx website for analysis.
write.table(sig_matrix,file = "RCC_signature_matrix.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##2 Making Mixture.
###Loading and processing TCGA-KIRC RNA-seq data.
setwd("E:/Raw Data/TCGA_data/TCGA-KIRC")
load("TCGA-KIRC_mRNA.Rdata") #Downloaded from the Genomic Data Commons (GDC) database.
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
t10 <- t10[, -c(1, 615)]
t10 <- as.matrix(t10)
saveRDS(t10, file = "TCGA-KIRC_RNAseq_TPM.rds")

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
saveRDS(t11, file = "TCGA-KIRC_tumor_RNAseq_TPM.rds")
saveRDS(t12, file = "TCGA-KIRC_normal_RNAseq_TPM.rds")

###The genes in the TCGA-KIRC tumor sample matrix aligned with the genes in the Seurat object matrix.
t11 <- as.data.frame(t11)
Gene <- rownames(RCC@assays$RNA@meta.features)
t13 <- t11[Gene,]
t13 <- t13[!rowSums(is.na(t13)),]

###Save the TCGA-KIRC Mixture.
saveRDS(t13, file = "TCGA-KIRC_tumor_RNAseq_TPM_filterd.rds")

###Save the TCGA-KIRC Mixture in txt format and upload them to the CIBERSORTx website for analysis.
write.table(t13,file = "TCGA-KIRC_Mixture.txt",quote = F,sep = "\t",row.names = T,col.names = T)

##3 Run CIBERSORTx. Detailed steps are in the CIBERSORTx.pdf file. #In the first row and first column of the 'txt format' Mixture, write 'Gene', separated from the rest of the first row by TAB. Similarly, in the 'txt format' signature matrix, write 'NAME' in the first row and first column, separated by TAB from the rest of the first row. Then, save the file again.

##4 Survival analysis.
###loading the CIBERSORTx result.
setwd("E:/Raw Data/cibersort_results/TCGA-KIRC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-KIRC_Results.csv")

###Loading TCGA-KIRC clinical data.
setwd("E:/Raw Data/TCGA_data/TCGA-KIRC")
clin2 <- readRDS(file = "TCGA-KIRC_clinical.rds") #Downloaded from the Genomic Data Commons (GDC) database.

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
t17 <- t17[, c(1, 7, 10, 11)]

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
saveRDS(t24, file = "TCGA-KIRC_clinical_processed_new.rds")

#Retain only the necessary clinical information.
t24 <- t24[, c(1:3)]

###CIBERSORTx results are integrated with clinical data.
t25 <- ciber[ciber$Mixture %in% t24$barcode,]
t26 <- merge(t24, t25, by.x = "barcode", by.y = "Mixture")
rownames(t26) <- t26[,1]
saveRDS(t26, file = "TCGA-KIRC_clinical_processed_cibersort_score.rds")

###Processing integrated data. Retain only clinical information and apCAFs signature scores.
t27 <- t26[, c(1, 2, 3, 4)]

###Exclude samples where OS.time is less than or equal to 0.
s1 <- t27
s2 <- s1[order(s1$OS.time),]
s3 <- s2[-c(1:11), ]

###Arrange in increasing order according to apCAFs scores.
s4 <- s3[order(s3$apCAFs),]

###Survival analysis.
library(survminer)
library(survival)

###Grouping.
s5<- s4[c(1:259),]
s5$group <- "low"
s6<- s4[c(260:518),]
s6$group <- "high"
s7 <- rbind(s5, s6)
table(s7$group)

###Survivorship curve.
sfit <- survfit(Surv(OS.time, status)~group, data=s7)

###Visualization.
#KM survival curve for Supplementary Figure 7G.
survp <- ggsurvplot(sfit,conf.int = T,pval.method = T, pval = T, legend.title='apCAFs signature score', legend.labs=c('High n = 259','Low  n = 259'), legend = c(0.75, 0.9), risk.table = F, xlab='days', ylab='Survival probability', palette = "jco", font.legend = list(size = 15))
pdf("ggsurvplot_TCGA-KIRC_apCAFs_signature_cibersortx.pdf",width = 5,height = 4)
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
shapiro.test(data$CD8..T.cells)

###Visualization.
#scatter plot for Supplementary Figure 7N.
pdf(file = "TCGA-KIRC_cibersortx_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
