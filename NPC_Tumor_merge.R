#NPC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of NPC tumors.

#Load Seurat objects with annotated cell types. 
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC <- readRDS(file = "NPC_resolution0.2_annotation_sub.rds") #The detailed code for NPC_resolution0.2_annotation_sub.rds is at line 616 of the Supplementary Figure 5A-B.R file.
NPC$celltype <- NPC$Celltype
NPC_CAF <- readRDS(file = "NPC_Fib_T_res0.1_annotation_sub.rds") #The detailed code for NPC_Fib_T_res0.1_annotation_sub.rds is at line 771 of the Supplementary Figure 5A-B.R file.
NPC_M <- readRDS(file = "NPC_Myeloid_res0.3_annotation_sub.rds") #The detailed code for NPC_Myeloid_res0.3_annotation_sub.rds is at line 187 of the Supplementary Figure 6E-F.R file.
NPC_T <- readRDS(file = "NPC_T_cells_res0.2_annotation_sub.rds") #The detailed code for NPC_T_cells_res0.2_annotation_sub.rds is at line 347 of the Supplementary Figure 6E-F.R file.

#Just tumor samples.
NPC_tumor <- subset(NPC, subset = Sample %in% c("Tumor"))
NPC_M_tumor <- subset(NPC_M, subset = Sample %in% c("Tumor"))
NPC_T_tumor <- subset(NPC_T, subset = Sample %in% c("Tumor"))

#Check cell types.
table(NPC_tumor$celltype)
table(NPC_M_tumor$celltype)
table(NPC_T_tumor$celltype)
table(NPC_CAF$celltype)

#Retain major cell types without further clustering annotation.
NPC_other_tumor <- subset(NPC_tumor, subset = celltype %in% c("B cells", "Plasma cells", "pDCs", "Mast cells", "Endothelial cells", "Epithelial cells"))

#Convert epithelial cells to tumor cells.
Idents(NPC_other_tumor) <- NPC_other_tumor$celltype
NPC_other_tumor <- RenameIdents(NPC_other_tumor, "Epithelial cells" = "Tumor cells")
table(Idents(NPC_other_tumor))
NPC_other_tumor$celltype <- Idents(NPC_other_tumor)
table(NPC_other_tumor$celltype)

#Reintegrate into a single-cell RNA-seq Seurat object annotated with various cell types for NPC tumor samples.
NPC_merge <- merge(x = NPC_other_tumor, y = c(NPC_M_tumor, NPC_T_tumor, NPC_CAF))
NPC_merge$seurat_clusters <- as.factor(as.character(NPC_merge$seurat_clusters))
NPC_merge$Celltype <- as.factor(as.character(NPC_merge$Celltype))
NPC_merge$orig.ident <- as.factor(as.character(NPC_merge$orig.ident))
NPC_merge$Sample <- as.factor(as.character(NPC_merge$Sample))
NPC_merge$celltype <- as.factor(as.character(NPC_merge$celltype))

#Check cell types and sample types.
table(NPC_merge$celltype)
table(NPC_merge$Sample)
saveRDS(NPC_merge, file = "NPC_Tumor_merge.rds")
