#CRC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of CRC tumors.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
CRC <- readRDS(file = "CRC_HRA000979_resolution0.1_annotation_sub.rds") #The detailed code for CRC_HRA000979_resolution0.1_annotation_sub.rds is at line 269 of the Supplementary Figure 5E-F.R file.
CRC$celltype <- CRC$Celltype
CRC_CAF <- readRDS(file = "CRC_HRA000979_Fib_T_res0.2_annotation.rds") #The detailed code for CRC_HRA000979_Fib_T_res0.2_annotation.rds is at line 411 of the Supplementary Figure 5E-F.R file.
CRC_M <- readRDS(file = "CRC_HRA000979_Myeloid_res0.2_annotation_sub.rds") #The detailed code for CRC_HRA000979_Myeloid_res0.2_annotation_sub.rds is at line 183 of the Supplementary Figure 6I-J.R file.
CRC_T <- readRDS(file = "CRC_HRA000979_T_cells_res0.1_annotation_sub.rds") #The detailed code for CRC_HRA000979_T_cells_res0.1_annotation_sub.rds is at line 345 of the Supplementary Figure 6I-J.R file.

#Just tumor samples.
CRC_tumor <- subset(CRC, subset = Sample %in% c("Tumor"))
CRC_M_tumor <- subset(CRC_M, subset = Sample %in% c("Tumor"))
CRC_T_tumor <- subset(CRC_T, subset = Sample %in% c("Tumor"))

#Check cell types.
table(CRC_tumor$celltype)
table(CRC_M_tumor$celltype)
table(CRC_T_tumor$celltype)
table(CRC_CAF$celltype)

#Retain major cell types without further clustering annotation.
CRC_other_tumor <- subset(CRC_tumor, subset = celltype %in% c("B cells", "Plasma cells", "Mast cells", "Endothelial cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(CRC_other_tumor) <- CRC_other_tumor$celltype
CRC_other_tumor <- RenameIdents(CRC_other_tumor, "Epithelial cells" = "Tumor cells")
table(Idents(CRC_other_tumor))
CRC_other_tumor$celltype <- Idents(CRC_other_tumor)
table(CRC_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for CRC tumor samples.
CRC_merge <- merge(x = CRC_other_tumor, y = c(CRC_M_tumor, CRC_T_tumor, CRC_CAF))
CRC_merge$seurat_clusters <- as.factor(as.character(CRC_merge$seurat_clusters))
CRC_merge$Celltype <- as.factor(as.character(CRC_merge$Celltype))
CRC_merge$orig.ident <- as.factor(as.character(CRC_merge$orig.ident))
CRC_merge$Sample <- as.factor(as.character(CRC_merge$Sample))
CRC_merge$celltype <- as.factor(as.character(CRC_merge$celltype))

#Check cell types and sample types.
table(CRC_merge$celltype)
table(CRC_merge$Sample)
saveRDS(CRC_merge, file = "CRC_Tumor_merge.rds")
