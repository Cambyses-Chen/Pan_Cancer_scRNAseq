###codes for Supplementary Figure 9H-N.

library(Seurat)
library(ggplot2)
library(scales)

#1 Supplementary Figure 9H.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib_new <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib_new$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9H.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="NPC_Fib_CTSS.pdf",width=8, height=5)
DotPlot(NPC_Fib_new, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 9I.
#OV.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in OV.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV_CAFs <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The codes for OV_Fib_res0.2_annotation_sub.rds is located at line 406 of Figure 1G-J.R.
table(OV_CAFs$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9I.
CTSS = CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="OV_Fib_CTSS.pdf",width=8, height=5)
DotPlot(OV_CAFs, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#3 Supplementary Figure 9J.
#CRC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in CRC.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
CRC_Fib_new <- readRDS(file = "CRC_Fib_new.rds") #The codes for CRC_Fib_new.rds is located at line 178 of Supplementary Figure 9A-E.R.
table(CRC_Fib_new$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9J.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="CRC_Fib_CTSS.pdf",width=8, height=5)
DotPlot(CRC_Fib_new, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#4 Supplementary Figure 9K.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_CAFs <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_CAFs$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9K.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="BRCA_Fib_CTSS.pdf",width=8, height=5)
DotPlot(BRCA_CAFs, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#5 Supplementary Figure 9L.
#AM.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in AM.
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
AM_CAFs <- readRDS(file = "AM_Fib_res0.1_annotation_sub.rds") #The codes for AM_Fib_res0.1_annotation_sub.rds is located at line 577 of Supplementary Figure 5I-J.R.
table(AM_CAFs$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9L.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="AM_Fib_CTSS.pdf",width=8, height=5)
DotPlot(AM_CAFs, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#6 Supplementary Figure 9M.
#CM.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in CM.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM_CAFs <- readRDS(file = "CM_Fib_res0.4_annotation_sub.rds") #The codes for CM_Fib_res0.4_annotation_sub.rds is located at line 392 of Supplementary Figure 5K-L.R.
table(CM_CAFs$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9M.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="CM_Fib_CTSS.pdf",width=8, height=5)
DotPlot(CM_CAFs, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#7 Supplementary Figure 9N.
#RCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in RCC.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC_Fib_new <- readRDS(file = "RCC_Fib_new.rds") #The codes for RCC_Fib_new.rds is located at line 159 of Supplementary Figure 9F-G.R.
table(RCC_Fib_new$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Supplementary Figure 9N.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="RCC_Fib_CTSS.pdf",width=8, height=5)
DotPlot(RCC_Fib_new, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
