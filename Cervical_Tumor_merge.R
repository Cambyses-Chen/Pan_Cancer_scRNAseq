#Cervical_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of cervical cancers.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
Cervical <- readRDS(file = "Cervical_resolution0.3_annotation_2.rds") #The detailed code for Cervical_resolution0.3_annotation_2.rds is at line 268 of the Supplementary Figure 5C-D.R file.
Cervical$celltype <- Cervical$Celltype
Cervical_CAF <- readRDS(file = "Cervical_Fib_T_res1_annotation_sub.rds") #The detailed code for Cervical_Fib_T_res1_annotation_sub.rds is at line 418 of the Supplementary Figure 5C-D.R file.
Cervical_M <- readRDS(file = "Cervical_Myeloid_res0.2_annotation_sub.rds") #The detailed code for Cervical_Myeloid_res0.2_annotation_sub.rds is at line 186 of the Supplementary Figure 6G-H.R file.
Cervical_T <- readRDS(file = "Cervical_T_cells_res0.1_annotation_sub.rds") #The detailed code for Cervical_T_cells_res0.1_annotation_sub.rds is at line 344 of the Supplementary Figure 6G-H.R file.

#Just tumor samples.
Cervical_tumor <- subset(Cervical, subset = Sample %in% c("Tumor"))
Cervical_M_tumor <- subset(Cervical_M, subset = Sample %in% c("Tumor"))
Cervical_T_tumor <- subset(Cervical_T, subset = Sample %in% c("Tumor"))

#Check cell types.
table(Cervical_tumor$celltype)
table(Cervical_M_tumor$celltype)
table(Cervical_T_tumor$celltype)
table(Cervical_CAF$celltype)

#Retain major cell types without further clustering annotation.
Cervical_other_tumor <- subset(Cervical_tumor, subset = celltype %in% c("B cells", "Plasma cells", "NK cells", "Mast cells", "Endothelial cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(Cervical_other_tumor) <- Cervical_other_tumor$celltype
Cervical_other_tumor <- RenameIdents(Cervical_other_tumor, "Epithelial cells" = "Tumor cells")
table(Idents(Cervical_other_tumor))
Cervical_other_tumor$celltype <- Idents(Cervical_other_tumor)
table(Cervical_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for Cervical tumor samples.
Cervical_merge <- merge(x = Cervical_other_tumor, y = c(Cervical_M_tumor, Cervical_T_tumor, Cervical_CAF))
Cervical_merge$seurat_clusters <- as.factor(as.character(Cervical_merge$seurat_clusters))
Cervical_merge$Celltype <- as.factor(as.character(Cervical_merge$Celltype))
Cervical_merge$orig.ident <- as.factor(as.character(Cervical_merge$orig.ident))
Cervical_merge$Sample <- as.factor(as.character(Cervical_merge$Sample))
Cervical_merge$celltype <- as.factor(as.character(Cervical_merge$celltype))

#Check cell types and sample types.
table(Cervical_merge$celltype)
table(Cervical_merge$Sample)
saveRDS(Cervical_merge, file = "Cervical_Tumor_merge.rds")
