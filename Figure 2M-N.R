###codes for Figure 2M-N.

library(Seurat)
library(ggplot2)
library(scales)

#1 Figure 2M.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib_new <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib_new$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Figure 2M.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="HNSCC_Fib_CTSS.pdf",width=8, height=5)
DotPlot(HNSCC_Fib_new, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Figure 2N.
#CC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in CC.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CC_Fib_new <- readRDS(file = "CC_Fib_new.rds") #The codes for CC_Fib_new.rds is located at line 118 of Supplementary Figure 9A-E.R.
table(CC_Fib_new$celltype)

###Visualization.
#The expression status of classical antigen processing and presentation-related molecular components for Figure 2N.
CTSS = c("TAP2", "TAP1", "RAB11A", "LNPEP", "CTSS")
pdf(file="CC_Fib_CTSS.pdf",width=8, height=5)
DotPlot(CC_Fib_new, features = CTSS, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
