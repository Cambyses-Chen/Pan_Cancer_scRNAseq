###codes for Supplementary Figure 9A-E.

library(Seurat)
library(ggplot2)
library(scales)

#1 Supplementary Figure 9A.
#NPC.
#Load single-cell RNA sequencing data of Fibroblasts major cell type in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib.rds")  #The codes for NPC_Fib.rds is located at line 645 of Supplementary Figure 5A-B.R.
table(NPC_Fib$Sample)
table(NPC_Fib$Celltype)

#Obtain single-cell RNA sequencing data of normal fibroblasts related to NPC.
NPC_NFs <- subset(NPC_Fib, subset = Sample %in% c("Normal"))
NPC_NFs$Sample <- as.factor(as.character(NPC_NFs$Sample))
NPC_NFs$orig.ident <- as.factor(as.character(NPC_NFs$orig.ident))
table(NPC_NFs$Sample)
table(NPC_NFs$Celltype)
NPC_NFs$celltype <- "NFs"
table(NPC_NFs$celltype)

#Load annotated single-cell RNA sequencing data of CAFs subpopulations in NPC.
NPC_CAFs <- readRDS(file = "NPC_Fib_T_res0.1_annotation_sub.rds")  #The codes for NPC_Fib_T_res0.1_annotation_sub.rds is located at line 771 of Supplementary Figure 5A-B.R.

#Merge into new annotated single-cell RNA sequencing data of fibroblast subpopulations in NPC.
NPC_Fib_new <- merge(x = NPC_CAFs, y = NPC_NFs)

#Ranking.
Idents(NPC_Fib_new) <- NPC_Fib_new$celltype
NPC_Fib_new <- RenameIdents(NPC_Fib_new, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs", "NFs" = "NFs")
table(Idents(NPC_Fib_new))
NPC_Fib_new$celltype <- Idents(NPC_Fib_new)
table(NPC_Fib_new$celltype)

###Normalizing the data.
NPC_Fib_new <- NormalizeData(NPC_Fib_new, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC_Fib_new <- FindVariableFeatures(NPC_Fib_new, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC_Fib_new)
NPC_Fib_new <- ScaleData(NPC_Fib_new, features = all.genes)
NPC_Fib_new <- ScaleData(NPC_Fib_new, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC_Fib_new <- RunPCA(NPC_Fib_new, features = VariableFeatures(object = NPC_Fib_new))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC_Fib_new <- RunHarmony(NPC_Fib_new, group.by.vars = "orig.ident")
NPC_Fib_new <- RunUMAP(NPC_Fib_new, reduction = "harmony", dims = 1:20)
NPC_Fib_new <- RunTSNE(NPC_Fib_new, reduction = "harmony", dims = 1:20)
names(NPC_Fib_new@reductions)
saveRDS(NPC_Fib_new, file = "NPC_Fib_new.rds")

###Visualization.
#The expression levels of C1Q molecules for Supplementary Figure 9A.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="NPC_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(NPC_Fib_new, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 9B.
#CC.
#Load single-cell RNA sequencing data of Fibroblasts major cell type in CC (Cervical Cancers). 
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CC_Fib <- readRDS(file = "Cervical_Fib.rds")  #The codes for Cervical_Fib.rds is located at line 296 of Supplementary Figure 5C-D.R.
table(CC_Fib$Sample)
table(CC_Fib$Celltype)

#Obtain single-cell RNA sequencing data of normal fibroblasts related to CC.
CC_NFs <- subset(CC_Fib, subset = Sample %in% c("Normal"))
CC_NFs$Sample <- as.factor(as.character(CC_NFs$Sample))
CC_NFs$orig.ident <- as.factor(as.character(CC_NFs$orig.ident))
table(CC_NFs$Sample)
table(CC_NFs$Celltype)
CC_NFs$celltype <- "NFs"
table(CC_NFs$celltype)

#Load annotated single-cell RNA sequencing data of CAFs subpopulations in CC.
CC_CAFs <- readRDS(file = "Cervical_Fib_T_res1_annotation_sub.rds")  #The codes for Cervical_Fib_T_res1_annotation_sub.rds is located at line 418 of Supplementary Figure 5C-D.R.

#Merge into new annotated single-cell RNA sequencing data of fibroblast subpopulations in CC.
CC_Fib_new <- merge(x = CC_CAFs, y = CC_NFs)

#Ranking.
Idents(CC_Fib_new) <- CC_Fib_new$celltype
CC_Fib_new <- RenameIdents(CC_Fib_new, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs", "NFs" = "NFs")
table(Idents(CC_Fib_new))
CC_Fib_new$celltype <- Idents(CC_Fib_new)
table(CC_Fib_new$celltype)

###Normalizing the data.
CC_Fib_new <- NormalizeData(CC_Fib_new, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CC_Fib_new <- FindVariableFeatures(CC_Fib_new, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CC_Fib_new)
CC_Fib_new <- ScaleData(CC_Fib_new, features = all.genes)
CC_Fib_new <- ScaleData(CC_Fib_new, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CC_Fib_new <- RunPCA(CC_Fib_new, features = VariableFeatures(object = CC_Fib_new))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CC_Fib_new <- RunHarmony(CC_Fib_new, group.by.vars = "orig.ident")
CC_Fib_new <- RunUMAP(CC_Fib_new, reduction = "harmony", dims = 1:20)
CC_Fib_new <- RunTSNE(CC_Fib_new, reduction = "harmony", dims = 1:20)
names(CC_Fib_new@reductions)
saveRDS(CC_Fib_new, file = "CC_Fib_new.rds")

###Visualization.
#The expression levels of C1Q molecules for Supplementary Figure 9B.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="CC_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(CC_Fib_new, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#3 Supplementary Figure 9C.
#CRC.
#Load single-cell RNA sequencing data of Fibroblasts major cell type in CRC.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
CRC_Fib <- readRDS(file = "CRC_HRA000979_Fib.rds")  #The codes for CRC_HRA000979_Fib.rds is located at line 296 of Supplementary Figure 5E-F.R.
table(CRC_Fib$Sample)
table(CRC_Fib$Celltype)

#Obtain single-cell RNA sequencing data of normal fibroblasts related to CRC.
CRC_NFs <- subset(CRC_Fib, subset = Sample %in% c("Normal"))
CRC_NFs$Sample <- as.factor(as.character(CRC_NFs$Sample))
CRC_NFs$orig.ident <- as.factor(as.character(CRC_NFs$orig.ident))
table(CRC_NFs$Sample)
table(CRC_NFs$Celltype)
CRC_NFs$celltype <- "NFs"
table(CRC_NFs$celltype)

#Load annotated single-cell RNA sequencing data of CAFs subpopulations in CRC.
CRC_CAFs <- readRDS(file = "CRC_HRA000979_Fib_T_res0.2_annotation.rds")  #The codes for CRC_HRA000979_Fib_T_res0.2_annotation.rds is located at line 411 of Supplementary Figure 5E-F.R.

#Merge into new annotated single-cell RNA sequencing data of fibroblast subpopulations in CRC.
CRC_Fib_new <- merge(x = CRC_CAFs, y = CRC_NFs)

#Ranking.
Idents(CRC_Fib_new) <- CRC_Fib_new$celltype
CRC_Fib_new <- RenameIdents(CRC_Fib_new, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "MYH11+ CAFs" = "MYH11+ CAFs", "iCAFs" = "iCAFs", "NFs" = "NFs")
table(Idents(CRC_Fib_new))
CRC_Fib_new$celltype <- Idents(CRC_Fib_new)
table(CRC_Fib_new$celltype)

###Normalizing the data.
CRC_Fib_new <- NormalizeData(CRC_Fib_new, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CRC_Fib_new <- FindVariableFeatures(CRC_Fib_new, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CRC_Fib_new)
CRC_Fib_new <- ScaleData(CRC_Fib_new, features = all.genes)
CRC_Fib_new <- ScaleData(CRC_Fib_new, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CRC_Fib_new <- RunPCA(CRC_Fib_new, features = VariableFeatures(object = CRC_Fib_new))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CRC_Fib_new <- RunHarmony(CRC_Fib_new, group.by.vars = "orig.ident")
CRC_Fib_new <- RunUMAP(CRC_Fib_new, reduction = "harmony", dims = 1:20)
CRC_Fib_new <- RunTSNE(CRC_Fib_new, reduction = "harmony", dims = 1:20)
names(CRC_Fib_new@reductions)
saveRDS(CRC_Fib_new, file = "CRC_Fib_new.rds")

###Visualization.
#The expression levels of C1Q molecules for Supplementary Figure 9C.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="CRC_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(CRC_Fib_new, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#4 Supplementary Figure 9D.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_CAFs <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_CAFs$celltype)

###Visualization.
#The expression levels of C1Q molecules for Supplementary Figure 9D.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="BRCA_Fib_Co_stimulation.pdf",width=6, height=5)
DotPlot(BRCA_CAFs, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#5 Supplementary Figure 9E.
#AM.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in AM.
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
AM_CAFs <- readRDS(file = "AM_Fib_res0.1_annotation_sub.rds") #The codes for AM_Fib_res0.1_annotation_sub.rds is located at line 577 of Supplementary Figure 5I-J.R.
table(AM_CAFs$celltype)

###Visualization.
#The expression levels of C1Q molecules for Supplementary Figure 9E.
Co_stimulation = c("C1QC", "C1QB", "C1QA")
pdf(file="AM_Fib_Co_stimulation.pdf",width=6, height=5.5)
DotPlot(AM_CAFs, features = Co_stimulation, dot.scale = 10)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
