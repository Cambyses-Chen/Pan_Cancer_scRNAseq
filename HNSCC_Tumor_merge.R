#HNSCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of HNSCC tumors.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC <- readRDS(file = "HNSCC_resolution0.3_annotation_sub.rds") #The detailed code for HNSCC_resolution0.3_annotation_sub.rds is at line 273 of the Figure 1C-F.R file.
HNSCC$celltype <- HNSCC$Celltype
HNSCC_CAF <- readRDS(file = "HNSCC_Fib_T_res0.2_annotation_sub.rds") #The detailed code for HNSCC_Fib_T_res0.2_annotation_sub.rds is at line 418 of the Figure 1C-F.R file.
HNSCC_M <- readRDS(file = "HNSCC_Myeloid_res0.4_annotation_sub.rds") #The detailed code for HNSCC_Myeloid_res0.4_annotation_sub.rds is at line 189 of the Supplementary Figure 6A-B.R file.
HNSCC_T <- readRDS(file = "HNSCC_T_cells_res0.1_annotation_sub.rds") #The detailed code for HNSCC_T_cells_res0.1_annotation_sub.rds is at line 347 of the Supplementary Figure 6A-B.R file.

#Just tumor samples.
HNSCC_tumor <- subset(HNSCC, subset = Sample %in% c("Primary", "Metastatic"))
HNSCC_M_tumor <- subset(HNSCC_M, subset = Sample %in% c("Primary", "Metastatic"))
HNSCC_T_tumor <- subset(HNSCC_T, subset = Sample %in% c("Primary", "Metastatic"))

#Check cell types.
table(HNSCC_tumor$celltype)
table(HNSCC_M_tumor$celltype)
table(HNSCC_T_tumor$celltype)
table(HNSCC_CAF$celltype)

#Retain major cell types without further clustering annotation.
HNSCC_other_tumor <- subset(HNSCC_tumor, subset = celltype %in% c("B cells", "Plasma cells", "pDCs", "Mast cells", "Endothelial cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(HNSCC_other_tumor) <- HNSCC_other_tumor$celltype
HNSCC_other_tumor <- RenameIdents(HNSCC_other_tumor, "Epithelial cells" = "Tumor cells")
table(Idents(HNSCC_other_tumor))
HNSCC_other_tumor$celltype <- Idents(HNSCC_other_tumor)
table(HNSCC_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for HNSCC tumor samples.
HNSCC_merge <- merge(x = HNSCC_other_tumor, y = c(HNSCC_M_tumor, HNSCC_T_tumor, HNSCC_CAF))
HNSCC_merge$seurat_clusters <- as.factor(as.character(HNSCC_merge$seurat_clusters))
HNSCC_merge$Celltype <- as.factor(as.character(HNSCC_merge$Celltype))
HNSCC_merge$orig.ident <- as.factor(as.character(HNSCC_merge$orig.ident))
HNSCC_merge$Sample <- as.factor(as.character(HNSCC_merge$Sample))
HNSCC_merge$celltype <- as.factor(as.character(HNSCC_merge$celltype))

#Check cell types and sample types.
table(HNSCC_merge$celltype)
table(HNSCC_merge$Sample)
saveRDS(HNSCC_merge, file = "HNSCC_Tumor_merge.rds")
