###codes for Figure 6G-L.

library(Seurat)
library(ggplot2)
library(scales)
library(ggpubr)
library(nortest)

#1 Figure 6G.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

###Extract the raw count matrix.
expr <- as.data.frame(HNSCC_Fib@assays$RNA@counts)

###Refine the raw count matrix by selecting genes using scFEA.human.genes.
setwd("E:/Raw Data/FLUXestimator_data")
gene <- read.delim(file = "scFEA.human.genes.txt") #scFEA.human.genes.txt was downloaded from the FLUXestimator website.
list <- as.list(gene)
gene2 <- list[[1]]
expr2 <- expr[gene2, ]
expr3 <- expr2[!rowSums(is.na(expr2)),]

###Save the refined raw count matrix in csv format and upload it to the FLUXestimator website for analysis.
setwd("E:/Raw Data/FLUXestimator_data/HNSCC")
write.csv(expr3, file = "HNSCC_Fib_expr_counts.csv",quote = F,row.names = T,col.names = T)

###Run scFEA in the FLUXestimator website. ###The detailed steps are in a PDF file called FLUXestimator.pdf.

###After running scFEA, read scFEA.Glucose-TCAcycle.human.moduleinfo and predicted metabolic flux result.
setwd("E:/Raw Data/FLUXestimator_data/HNSCC")
t <- read.csv(file = "scFEA.Glucose-TCAcycle.human.moduleinfo.csv")###Get module Information. #scFEA.Glucose-TCAcycle.human.moduleinfo.csv was downloaded from the FLUXestimator website.
data <- read.csv('HNSCC_Fib_expr_counts_flux.csv')

###Extract cell types information from the Seurat object.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)
ID <- HNSCC_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###predicted metabolic flux result are integrated with cell types information.
data2 <- merge(ID2, data, by.x = "barcode", by.y = "X")
rownames(data2) <- data2[, 1]

###Anderson-Darling test for normality.
ad.test(data2$M_6) #M_6 represent Pyruvate -> Lactate. The information is sourced from scFEA.Glucose-TCAcycle.human.moduleinfo.csv.

###Visualization.
setwd("E:/Raw Data/FLUXestimator_data/HNSCC")
table(data2$celltype)

#boxplot for Figure 6G.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"), c("apCAFs", "NFs") )
pdf(file = "ggboxplot_HNSCC_Pyruvate_Lactate.pdf", width=5, height=4.5)
ggboxplot(data2, x = "celltype", y = "M_6", color = "celltype", 
          add = "jitter", legend = "none", title = "Pyruvate -> Lactate", xlab = "", ylab = "Estimated Flux Level", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.65, 0.7, 0.75, 0.8))
dev.off()

#2 Figure 6H.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

###Extract the raw count matrix.
expr <- as.data.frame(NPC_Fib@assays$RNA@counts)

###Refine the raw count matrix by selecting genes using scFEA.human.genes.
setwd("E:/Raw Data/FLUXestimator_data")
gene <- read.delim(file = "scFEA.human.genes.txt") #scFEA.human.genes.txt was downloaded from the FLUXestimator website.
list <- as.list(gene)
gene2 <- list[[1]]
expr2 <- expr[gene2, ]
expr3 <- expr2[!rowSums(is.na(expr2)),]

###Save the refined raw count matrix in csv format and upload it to the FLUXestimator website for analysis.
setwd("E:/Raw Data/FLUXestimator_data/NPC")
write.csv(expr3, file = "NPC_Fib_expr_counts.csv",quote = F,row.names = T,col.names = T)

###Run scFEA in the FLUXestimator website. ###The detailed steps are in a PDF file called FLUXestimator.pdf.

###After running scFEA, read scFEA.Glucose-TCAcycle.human.moduleinfo and predicted metabolic flux result.
setwd("E:/Raw Data/FLUXestimator_data/NPC")
t <- read.csv(file = "scFEA.Glucose-TCAcycle.human.moduleinfo.csv")###Get module Information. #scFEA.Glucose-TCAcycle.human.moduleinfo.csv was downloaded from the FLUXestimator website.
data <- read.csv('NPC_Fib_expr_counts_flux.csv')

###Extract cell types information from the Seurat object.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)
ID <- NPC_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###predicted metabolic flux results are integrated with cell types information.
data2 <- merge(ID2, data, by.x = "barcode", by.y = "X")
rownames(data2) <- data2[, 1]

###Shapiro-Wilk test for normality.
shapiro.test(data2$M_6)#M_6 represent Pyruvate -> Lactate. The information is sourced from scFEA.Glucose-TCAcycle.human.moduleinfo.csv.

###Visualization.
setwd("E:/Raw Data/FLUXestimator_data/NPC")
table(data2$celltype)

#boxplot for Figure 6H.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"), c("apCAFs", "NFs") )
pdf(file = "ggboxplot_NPC_Pyruvate_Lactate.pdf", width=5, height=4.5)
ggboxplot(data2, x = "celltype", y = "M_6", color = "celltype", 
          add = "jitter", legend = "none", title = "Pyruvate -> Lactate", xlab = "", ylab = "Estimated Flux Level", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.65, 0.7, 0.75, 0.8))
dev.off()

#3 Figure 6I.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_Fib$celltype)

###Extract the raw count matrix.
expr <- as.data.frame(BRCA_Fib@assays$RNA@counts)

###Refine the raw count matrix by selecting genes using scFEA.human.genes.
setwd("E:/Raw Data/FLUXestimator_data")
gene <- read.delim(file = "scFEA.human.genes.txt") #scFEA.human.genes.txt was downloaded from the FLUXestimator website.
list <- as.list(gene)
gene2 <- list[[1]]
expr2 <- expr[gene2, ]
expr3 <- expr2[!rowSums(is.na(expr2)),]

###Save the refined raw count matrix in csv format and upload it to the FLUXestimator website for analysis.
setwd("E:/Raw Data/FLUXestimator_data/BRCA")
write.csv(expr3, file = "BRCA_Fib_expr_counts.csv",quote = F,row.names = T,col.names = T)

###Run scFEA in the FLUXestimator website. ###The detailed steps are in a PDF file called FLUXestimator.pdf.

###After running scFEA, read scFEA.Glucose-TCAcycle.human.moduleinfo and predicted metabolic flux result.
setwd("E:/Raw Data/FLUXestimator_data/BRCA")
t <- read.csv(file = "scFEA.Glucose-TCAcycle.human.moduleinfo.csv")###Get module Information. #scFEA.Glucose-TCAcycle.human.moduleinfo.csv was downloaded from the FLUXestimator website.
data <- read.csv('BRCA_Fib_expr_counts_flux.csv')

###Extract cell types information from the Seurat object.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_Fib$celltype)
ID <- BRCA_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###predicted metabolic flux result are integrated with cell types information.
data2 <- merge(ID2, data, by.x = "barcode", by.y = "X")
rownames(data2) <- data2[, 1]

###Anderson-Darling test for normality.
ad.test(data2$M_6)#M_6 represent Pyruvate -> Lactate. The information is sourced from scFEA.Glucose-TCAcycle.human.moduleinfo.csv.

###Visualization.
setwd("E:/Raw Data/FLUXestimator_data/BRCA")
table(data2$celltype)

#boxplot for Figure 6I.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"))
pdf(file = "ggboxplot_BRCA_Pyruvate_Lactate.pdf", width=4.2, height=4.5)
ggboxplot(data2, x = "celltype", y = "M_6", color = "celltype", 
          add = "jitter", legend = "none", title = "Pyruvate -> Lactate", xlab = "", ylab = "Estimated Flux Level", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.75, 0.85, 0.95, 1.05))
dev.off()

#4 Figure 6J.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

###Visualization.
#The expression levels of GLUT1 genes SLC2A1 for Figure 6J.
SLC2A1 = c("SLC2A1")
pdf(file="HNSCC_Fib_SLC2A1.pdf",width=7, height=5)
DotPlot(HNSCC_Fib, features = SLC2A1, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#5 Figure 6K.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

###Visualization.
#The expression levels of GLUT1 genes SLC2A1 for Figure 6K.
SLC2A1 = c("SLC2A1")
pdf(file="NPC_Fib_SLC2A1.pdf",width=7, height=5)
DotPlot(NPC_Fib, features = SLC2A1, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#6 Figure 6L.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_Fib$celltype)

###Visualization.
#The expression levels of GLUT1 genes SLC2A1 for Figure 6L.
SLC2A1 = c("SLC2A1")
pdf(file="BRCA_Fib_SLC2A1.pdf",width=7, height=5)
DotPlot(BRCA_Fib, features = SLC2A1, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
