#RCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of RCC tumors.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC <- readRDS(file = "RCC_resolution0.1_annotation2.rds") #The detailed code for RCC_resolution0.1_annotation2.rds is at line 288 of the Supplementary Figure 5M-N.R file.
RCC$celltype <- RCC$Celltype
RCC_CAF <- readRDS(file = "RCC_Fib_T_res0.3_annotation_sub.rds") #The detailed code for RCC_Fib_T_res0.3_annotation_sub.rds is at line 436 of the Supplementary Figure 5M-N.R file.
RCC_M <- readRDS(file = "RCC_Myeloid_res0.1_annotation2.rds") #The detailed code for RCC_Myeloid_res0.1_annotation2.rds is at line 176 of the Supplementary Figure 6Q-R.R file.
RCC_T <- readRDS(file = "RCC_T_cells_res0.1_annotation_sub.rds") #The detailed code for RCC_T_cells_res0.1_annotation_sub.rds is at line 334 of the Supplementary Figure 6Q-R.R file.

#Just tumor samples.
RCC_tumor <- subset(RCC, subset = Sample %in% c("Primary tumor"))
RCC_M_tumor <- subset(RCC_M, subset = Sample %in% c("Primary tumor"))
RCC_T_tumor <- subset(RCC_T, subset = Sample %in% c("Primary tumor"))

#Check cell types.
table(RCC_tumor$celltype)
table(RCC_M_tumor$celltype)
table(RCC_T_tumor$celltype)
table(RCC_CAF$celltype)

#Retain major cell types without further clustering annotation.
RCC_other_tumor <- subset(RCC_tumor, subset = celltype %in% c("B cells", "NK cells", "Mast cells", "Endothelial cells", "Clear cells"))

#Convert clear cells to tumor cells.
Idents(RCC_other_tumor) <- RCC_other_tumor$celltype
RCC_other_tumor <- RenameIdents(RCC_other_tumor, "Clear cells" = "Tumor cells")
table(Idents(RCC_other_tumor))
RCC_other_tumor$celltype <- Idents(RCC_other_tumor)
table(RCC_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for RCC tumor samples.
RCC_merge <- merge(x = RCC_other_tumor, y = c(RCC_M_tumor, RCC_T_tumor, RCC_CAF))
RCC_merge$seurat_clusters <- as.factor(as.character(RCC_merge$seurat_clusters))
RCC_merge$Celltype <- as.factor(as.character(RCC_merge$Celltype))
RCC_merge$orig.ident <- as.factor(as.character(RCC_merge$orig.ident))
RCC_merge$Sample <- as.factor(as.character(RCC_merge$Sample))
RCC_merge$celltype <- as.factor(as.character(RCC_merge$celltype))

#Check cell types and sample types.
table(RCC_merge$celltype)
table(RCC_merge$Sample)
saveRDS(RCC_merge, file = "RCC_Tumor_merge.rds")
