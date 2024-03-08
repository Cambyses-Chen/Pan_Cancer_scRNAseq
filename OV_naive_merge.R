#OV_naive_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of OV tumors without drug treatment.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV <- readRDS(file = "OV_resolution0.2_annotation2.rds") #The detailed code for OV_resolution0.2_annotation2.rds is at line 261 of the Figure 1G-J.R file.
OV$celltype <- OV$Celltype
OV_CAF <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The detailed code for OV_Fib_res0.2_annotation_sub.rds is at line 406 of the Figure 1G-J.R file.
OV_M <- readRDS(file = "OV_Myeloid_res0.4_annotation_sub.rds") #The detailed code for OV_Myeloid_res0.4_annotation_sub.rds is at line 187 of the Supplementary Figure 6C-D.R file.
OV_T <- readRDS(file = "OV_T_cells_res0.1_annotation_sub.rds") #The detailed code for OV_T_cells_res0.1_annotation_sub.rds is at line 344 of the Supplementary Figure 6C-D.R file.

#Just OV tumors samples without drug treatment.
OV_naive <- subset(OV, subset = Treatment %in% c("naive"))
OV_CAF_naive <- subset(OV_CAF, subset = Treatment %in% c("naive"))
OV_M_naive <- subset(OV_M, subset = Treatment %in% c("naive"))
OV_T_naive <- subset(OV_T, subset = Treatment %in% c("naive"))

#Check cell types.
table(OV_naive$celltype)
table(OV_M_naive$celltype)
table(OV_T_naive$celltype)
table(OV_CAF_naive$celltype)

#Retain major cell types without further clustering annotation.
OV_other_naive <- subset(OV_naive, subset = celltype %in% c("B cells", "Plasma cells", "Mast cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(OV_other_naive) <- OV_other_naive$celltype
OV_other_naive <- RenameIdents(OV_other_naive, "Epithelial cells" = "Tumor cells")
table(Idents(OV_other_naive))
OV_other_naive$celltype <- Idents(OV_other_naive)
table(OV_other_naive$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for OV tumors samples without drug treatment.
OV_merge <- merge(x = OV_other_naive, y = c(OV_M_naive, OV_T_naive, OV_CAF_naive))
OV_merge$seurat_clusters <- as.factor(as.character(OV_merge$seurat_clusters))
OV_merge$Celltype <- as.factor(as.character(OV_merge$Celltype))
OV_merge$orig.ident <- as.factor(as.character(OV_merge$orig.ident))
OV_merge$Sample <- as.factor(as.character(OV_merge$Sample))
OV_merge$Treatment <- as.factor(as.character(OV_merge$Treatment))
OV_merge$celltype <- as.factor(as.character(OV_merge$celltype))

#Check cell types and sample types.
table(OV_merge$celltype)
table(OV_merge$Treatment)
saveRDS(OV_merge, file = "OV_naive_merge.rds")
