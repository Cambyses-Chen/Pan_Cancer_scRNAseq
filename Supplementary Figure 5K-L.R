###codes for Supplementary Figure 5K-L.

#1 Supplementary Figure 5K.
#CM.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
CM1 <- Read10X_h5("GSM6622299_CM1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
CM2 <- Read10X_h5("GSM6622300_CM2_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
CM3 <- Read10X_h5("GSM6622301_CM3_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
CM1_lym <- Read10X_h5("GSM6622302_CM1_lym_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cm1 <- CreateSeuratObject(CM1, project = "CM1")
cm2 <- CreateSeuratObject(CM2, project = "CM2")
cm3 <- CreateSeuratObject(CM3, project = "CM3")
cm1_lym <- CreateSeuratObject(CM1_lym, project = "CM1_lym")

###Annotation.
cm1$Sample = "primary"
cm2$Sample = "primary"
cm3$Sample = "primary"
cm1_lym$Sample = "metastatic"

###Create one Seurat object.
CM <- merge(cm1, c(cm2, cm3, cm1_lym), add.cell.ids = c("CM1", "CM2", "CM3", "CM1_lym"))
CM
head(colnames(CM))
tail(colnames(CM))
table(CM$orig.ident)
saveRDS(CM, file = "CM_seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(CM, expression = nFeature_RNA > 200)
selected_f <- rownames(CM)[Matrix::rowSums(CM) > 3]
CM <- subset(CM, features = selected_f, cells = selected_c)
dim(CM)

###QC metrics and filtering.
CM[["percent.mt"]] <- PercentageFeatureSet(CM, pattern = "^MT-")
Idents(CM) <- CM$orig.ident
pdf(file="CM.featureViolin.pdf",width=6,height=6)           
VlnPlot(object = CM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="CM.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(CM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
CM <- subset(CM, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(CM)

###Visualization.
#VlnPlot for Supplementary Figure 1B.
Idents(CM) <- CM$orig.ident
pdf(file="CM.featureViolin_QC.pdf",width=5,height=6)           
VlnPlot(object = CM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
CM <- NormalizeData(CM, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CM <- FindVariableFeatures(CM, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CM)
CM <- ScaleData(CM, features = all.genes)
CM <- ScaleData(CM, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CM <- RunPCA(CM, features = VariableFeatures(object = CM))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CM <- RunHarmony(CM, group.by.vars = "orig.ident")
CM <- RunUMAP(CM, reduction = "harmony", dims = 1:20)
CM <- RunTSNE(CM, reduction = "harmony", dims = 1:20)
names(CM@reductions)

###Cluster the cells.
CM <- FindNeighbors(CM, reduction = "harmony", dims = 1:20)
saveRDS(CM, file = "CM_before_resolution.rds")

#Determine resolution.
test <- CM
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3H.
library(clustree)
library(patchwork)
pdf(file = "clustree_CM.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
CM <- FindClusters(CM, resolution = 0.2)
saveRDS(CM, file = "CM_resolution_0.2.rds")

###classic celltype marker expression across different clusters.
##1
celltype_marker=c(
  "EPCAM",# epithelial
  "PECAM1",#endothelial
  "COL3A1",# fibroblasts
  "CD163","AIF1",# myeloid
  "CD79A",#B cells
  "JCHAIN",# plasma cell
  "CD3D","CD8A","CD4",#T cells
  "GNLY","NKG7",#NK cells
  "PTPRC"#immune
)
pdf(file="CM_marker.pdf",width=44, height=33)
VlnPlot(CM, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="CM_NK.pdf",width=44, height=14)
VlnPlot(CM, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="CM_Epithelial.pdf",width=44, height=14)
VlnPlot(CM, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="CM_Plasma_cells.pdf",width=44, height=14)
VlnPlot(CM, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="CM_B_marker.pdf",width=44, height=19)
VlnPlot(CM, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="CM_Fibroblasts.pdf",width=44, height=14)
VlnPlot(CM, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="CM_Endothelial.pdf",width=44, height=5)
VlnPlot(CM, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="CM_T_marker.pdf",width=44, height=14)
VlnPlot(CM, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="CM_Mast.pdf",width=44, height=5)
VlnPlot(CM, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="CM_myeloid.pdf",width=44, height=9.5)
VlnPlot(CM, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="CM_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(CM, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="CM_pDC.pdf",width=44, height=20)
VlnPlot(CM, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="CM_neutrophil.pdf",width=44, height=5)
VlnPlot(CM, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="CM_melanocyte.pdf",width=44, height=30)
VlnPlot(CM, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="CM_clear_cell.pdf",width=44, height=10)
VlnPlot(CM, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(7, 8)
Myeloid_cells = c(9)
B_cells = c(4)
T_cells = c(1, 2)
Melanocytes = c(0, 3, 5)
Endothelial_cells = c(6)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         B_cells,
                         T_cells,
                         Melanocytes,
                         Endothelial_cells)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Melanocytes",length(Melanocytes)),
                     rep("Endothelial_cells",length(Endothelial_cells))
)

CM@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(CM@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CM@meta.data$Celltype)
table(CM$seurat_clusters)

###Rename.
Idents(CM) <- CM$Celltype
CM <- RenameIdents(CM, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "B_cells" = "B cells", "Myeloid_cells" = "Myeloid cells")
table(Idents(CM))
CM$Celltype <- Idents(CM)
table(CM$Celltype)
saveRDS(CM, file = "CM_resolution0.2_annotation.rds")

###Ranking.
Idents(CM) <- CM$Celltype
CM <- RenameIdents(CM, "T cells" = "T cells", "B cells" = "B cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Melanocytes" = "Melanocytes")
table(Idents(CM))
CM$Celltype <- Idents(CM)
table(CM$Celltype)
saveRDS(CM, file = "CM_resolution0.2_annotation2.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5K: left.
pdf(file="CM_umap_Celltype.pdf",width=5.8, height=4)
DimPlot(CM, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5K: right.
Idents(CM) <- 'Celltype'
My_levels <- c("T cells", "B cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Melanocytes")
Idents(CM) <- factor(Idents(CM), levels= My_levels)
features = c("MITF", "MLANA", "PMEL", "TYR", "VWF", "PECAM1", "COL3A1", "COL1A1", "LYZ", "AIF1", "MS4A1", "CD79A", "CD3D", "CD3E")
pdf(file="CM_features.pdf",width=8, height=7.5)
DotPlot(CM, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2H.
pdf(file="CM_umap_orig.ident.pdf",width=5.6, height=4)
DimPlot(CM, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5L.
###Cluster the Fibroblasts in tumor.
CM_Fib <- subset(CM, subset = Celltype %in% c("Fibroblasts"))
CM_Fib$Celltype <- as.factor(as.character(CM_Fib$Celltype))
CM_Fib$orig.ident <- as.factor(as.character(CM_Fib$orig.ident))
saveRDS(CM_Fib, file = "CM_Fib.rds")
Fib <- CM_Fib
table(Fib$Celltype)
table(Fib$Sample)

###Normalizing the data.
Fib <- NormalizeData(Fib, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
Fib <- FindVariableFeatures(Fib, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(Fib)
Fib <- ScaleData(Fib, features = all.genes)
Fib <- ScaleData(Fib, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
Fib <- RunPCA(Fib, features = VariableFeatures(object = Fib))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
Fib <- RunHarmony(Fib, group.by.vars = "orig.ident")
Fib <- RunUMAP(Fib, reduction = "harmony", dims = 1:20)
Fib <- RunTSNE(Fib, reduction = "harmony", dims = 1:20)
names(Fib@reductions)

###Cluster the cells.
Fib <- FindNeighbors(Fib, reduction = "harmony", dims = 1:20)
saveRDS(Fib, file = "CM_Fib_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4H.
library(clustree)
library(patchwork)
pdf(file = "clustree_CM_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.4.
Fib <- FindClusters(Fib, resolution = 0.4)
saveRDS(Fib, file = "CM_Fib_res0.4.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="CM_Fib_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="CM_Fib_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="CM_Fib_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="CM_Fib_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5", "MYH11")
pdf(file="CM_Fib_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(4)
eCAFs = c(5)
myCAFs = c(2)
iCAFs = c(1, 6)
unknown_CAFs = c(0, 3, 7, 8)
current.cluster.ids <- c(apCAFs,
                         eCAFs,
                         myCAFs,
                         iCAFs,
                         unknown_CAFs)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
                     rep("eCAFs",length(eCAFs)),
                     rep("myCAFs",length(myCAFs)),
                     rep("iCAFs",length(iCAFs)),
                     rep("unknown_CAFs",length(unknown_CAFs))
)

Fib@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Fib@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Fib$celltype)
table(Fib$seurat_clusters)

###Rename.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "unknown_CAFs" = "unknown CAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "CM_Fib_res0.4_annotation.rds")

###Exclude unknown cell subpopulations.
Fib <- subset(Fib, subset = celltype %in% c("apCAFs", "eCAFs", "myCAFs", "iCAFs"))
Fib$celltype <- as.factor(as.character(Fib$celltype))
Fib$orig.ident <- as.factor(as.character(Fib$orig.ident))
Fib$Sample <- as.factor(as.character(Fib$Sample))
table(Fib$celltype)

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "CM_Fib_res0.4_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5L: left.
pdf(file="CM_Fib_umap_celltype.pdf",width=5.2, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5L: right.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "eCAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "CXCL12", "CXCL1", "CCL2", "POSTN", "LOXL2", "MMP14", "MYH11", "RGS5", "ACTA2", "HLA-DRB1", "HLA-DPB1", "CD74")
pdf(file="CM_Fib_F_features_res0.4_annotation.pdf",width=8, height=6.5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
