#CM_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of CM tumors. 
#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM <- readRDS(file = "CM_resolution0.2_annotation2.rds") #The detailed code for CM_resolution0.2_annotation2.rds is at line 246 of the Supplementary Figure 5K-L.R file.
CM$celltype <- CM$Celltype
CM_CAF <- readRDS(file = "CM_Fib_res0.4_annotation_sub.rds") #The detailed code for CM_Fib_res0.4_annotation_sub.rds is at line 392 of the Supplementary Figure 5K-L.R file.
CM_M <- readRDS(file = "CM_Myeloid_res0.4_annotation_sub.rds") #The detailed code for CM_Myeloid_res0.4_annotation_sub.rds is at line 183 of the Supplementary Figure 6O-P.R file.
CM_T <- readRDS(file = "CM_T_cells_res0.2_annotation_sub.rds") #The detailed code for CM_T_cells_res0.2_annotation_sub.rds is at line 341 of the Supplementary Figure 6O-P.R file.

#Just tumor samples.
CM_tumor <- subset(CM, subset = Sample %in% c("primary", "metastatic"))
CM_M_tumor <- subset(CM_M, subset = Sample %in% c("primary", "metastatic"))
CM_T_tumor <- subset(CM_T, subset = Sample %in% c("primary", "metastatic"))

#Check cell types.
table(CM_tumor$celltype)
table(CM_M_tumor$celltype)
table(CM_T_tumor$celltype)
table(CM_CAF$celltype)

#Retain major cell types without further clustering annotation.
CM_other_tumor <- subset(CM_tumor, subset = celltype %in% c("B cells", "Endothelial cells", "Melanocytes"))

#Convert melanocytes to tumor cells.
Idents(CM_other_tumor) <- CM_other_tumor$celltype
CM_other_tumor <- RenameIdents(CM_other_tumor, "Melanocytes" = "Tumor cells")
table(Idents(CM_other_tumor))
CM_other_tumor$celltype <- Idents(CM_other_tumor)
table(CM_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for CM tumor samples.
CM_merge <- merge(x = CM_other_tumor, y = c(CM_M_tumor, CM_T_tumor, CM_CAF))
CM_merge$seurat_clusters <- as.factor(as.character(CM_merge$seurat_clusters))
CM_merge$Celltype <- as.factor(as.character(CM_merge$Celltype))
CM_merge$orig.ident <- as.factor(as.character(CM_merge$orig.ident))
CM_merge$Sample <- as.factor(as.character(CM_merge$Sample))
CM_merge$celltype <- as.factor(as.character(CM_merge$celltype))

#Check cell types and sample types.
table(CM_merge$celltype)
table(CM_merge$Sample)
saveRDS(CM_merge, file = "CM_Tumor_merge.rds")
