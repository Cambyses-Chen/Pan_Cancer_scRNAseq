#BC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of BRCA tumors.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BC <- readRDS(file = "BC_GSE176078_resolution0.1_annotation_sub.rds") #The detailed code for BC_GSE176078_resolution0.1_annotation_sub.rds is at line 308 of the Supplementary Figure 5G-H.R file.
BC$celltype <- BC$Celltype
BC_CAF <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The detailed code for BC_GSE176078_Fib_res0.2_annotation2.rds is at line 437 of the Supplementary Figure 5G-H.R file.
BC_M <- readRDS(file = "BC_GSE176078_Myeloid_res0.2_annotation_sub.rds") #The detailed code for BC_GSE176078_Myeloid_res0.2_annotation_sub.rds is at line 185 of the Supplementary Figure 6K-L.R file.
BC_T <- readRDS(file = "BC_GSE176078_T_cells_res0.1_annotation_sub.rds") #The detailed code for BC_GSE176078_T_cells_res0.1_annotation_sub.rds is at line 345 of the Supplementary Figure 6K-L.R file.

#Just tumor samples.
BC_tumor <- subset(BC, subset = Sample %in% c("ER", "HER2", "TNBC"))
BC_M_tumor <- subset(BC_M, subset = Sample %in% c("ER", "HER2", "TNBC"))
BC_T_tumor <- subset(BC_T, subset = Sample %in% c("ER", "HER2", "TNBC"))

#Check cell types.
table(BC_tumor$celltype)
table(BC_M_tumor$celltype)
table(BC_T_tumor$celltype)
table(BC_CAF$celltype)

#Retain major cell types without further clustering annotation.
BC_other_tumor <- subset(BC_tumor, subset = celltype %in% c("B cells", "Plasma cells", "Endothelial cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(BC_other_tumor) <- BC_other_tumor$celltype
BC_other_tumor <- RenameIdents(BC_other_tumor, "Epithelial cells" = "Tumor cells")
table(Idents(BC_other_tumor))
BC_other_tumor$celltype <- Idents(BC_other_tumor)
table(BC_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for BRCA tumor samples.
BC_merge <- merge(x = BC_other_tumor, y = c(BC_M_tumor, BC_T_tumor, BC_CAF))
BC_merge$seurat_clusters <- as.factor(as.character(BC_merge$seurat_clusters))
BC_merge$Celltype <- as.factor(as.character(BC_merge$Celltype))
BC_merge$orig.ident <- as.factor(as.character(BC_merge$orig.ident))
BC_merge$Sample <- as.factor(as.character(BC_merge$Sample))
BC_merge$celltype <- as.factor(as.character(BC_merge$celltype))

#Check cell types and sample types.
table(BC_merge$celltype)
table(BC_merge$Sample)
saveRDS(BC_merge, file = "BC_Tumor_merge.rds")
