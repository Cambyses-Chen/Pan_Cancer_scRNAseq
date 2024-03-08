###codes for Figure 2I-J.

library(Seurat)
library(ggplot2)
library(scales)

#1 Figure 2I.
#HNSCC.
#Load single-cell RNA sequencing data of Fibroblasts major cell type in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib.rds")  #The codes for HNSCC_Fib.rds is located at line 301 of Figure 1C-F.R.
table(HNSCC_Fib$Sample)
table(HNSCC_Fib$Celltype)

#Obtain single-cell RNA sequencing data of normal fibroblasts related to HNSCC.
HNSCC_NFs <- subset(HNSCC_Fib, subset = Sample %in% c("Normal"))
HNSCC_NFs$Sample <- as.factor(as.character(HNSCC_NFs$Sample))
HNSCC_NFs$orig.ident <- as.factor(as.character(HNSCC_NFs$orig.ident))
table(HNSCC_NFs$Sample)
table(HNSCC_NFs$Celltype)
HNSCC_NFs$celltype <- "NFs"
table(HNSCC_NFs$celltype)

#Load annotated single-cell RNA sequencing data of CAFs subpopulations in HNSCC.
HNSCC_CAFs <- readRDS(file = "HNSCC_Fib_T_res0.2_annotation_sub.rds")  #The codes for HNSCC_Fib_T_res0.2_annotation_sub.rds is located at line 418 of Figure 1C-F.R.

#Merge into new annotated single-cell RNA sequencing data of fibroblast subpopulations in HNSCC.
HNSCC_Fib_new <- merge(x = HNSCC_CAFs, y = HNSCC_NFs)

#Ranking.
Idents(HNSCC_Fib_new) <- HNSCC_Fib_new$celltype
HNSCC_Fib_new <- RenameIdents(HNSCC_Fib_new, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs", "NFs" = "NFs")
table(Idents(HNSCC_Fib_new))
HNSCC_Fib_new$celltype <- Idents(HNSCC_Fib_new)
table(HNSCC_Fib_new$celltype)

###Normalizing the data.
HNSCC_Fib_new <- NormalizeData(HNSCC_Fib_new, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
HNSCC_Fib_new <- FindVariableFeatures(HNSCC_Fib_new, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(HNSCC_Fib_new)
HNSCC_Fib_new <- ScaleData(HNSCC_Fib_new, features = all.genes)
HNSCC_Fib_new <- ScaleData(HNSCC_Fib_new, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
HNSCC_Fib_new <- RunPCA(HNSCC_Fib_new, features = VariableFeatures(object = HNSCC_Fib_new))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
HNSCC_Fib_new <- RunHarmony(HNSCC_Fib_new, group.by.vars = "orig.ident")
HNSCC_Fib_new <- RunUMAP(HNSCC_Fib_new, reduction = "harmony", dims = 1:20)
HNSCC_Fib_new <- RunTSNE(HNSCC_Fib_new, reduction = "harmony", dims = 1:20)
names(HNSCC_Fib_new@reductions)
saveRDS(HNSCC_Fib_new, file = "HNSCC_Fib_new.rds")

###Visualization.
#The expression levels of C1Q molecules for Figure 2I.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="HNSCC_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(HNSCC_Fib_new, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Figure 2J.
#OV.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in OV.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV_CAFs <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The codes for OV_Fib_res0.2_annotation_sub.rds is located at line 406 of Figure 1G-J.R.
table(OV_CAFs$celltype)

###Visualization.
#The expression levels of C1Q molecules for Figure 2J.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="OV_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(OV_CAFs, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
