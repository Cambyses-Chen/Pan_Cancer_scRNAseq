###codes for Supplementary Figure 6M-N.

#AM.
#1 Supplementary Figure 6M.
###Cluster the Myeloid cells.
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
AM <- readRDS(file = "AM_resolution0.1_annotation_sub.rds") #The detailed code for AM_resolution0.1_annotation_sub.rds is at line 436 of the Supplementary Figure 5I-J.R file.
AM_M <- subset(AM, subset = Celltype %in% c("Myeloid cells"))
AM_M$seurat_clusters <- as.factor(as.character(AM_M$seurat_clusters))
AM_M$Celltype <- as.factor(as.character(AM_M$Celltype))
AM_M$orig.ident <- as.factor(as.character(AM_M$orig.ident))
AM_M$Sample <- as.factor(as.character(AM_M$Sample))
table(AM_M$Celltype)
saveRDS(AM_M, file = "AM_Myeloid.rds")

###Normalizing the data.
AM_M <- NormalizeData(AM_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
AM_M <- FindVariableFeatures(AM_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(AM_M)
AM_M <- ScaleData(AM_M, features = all.genes)
AM_M <- ScaleData(AM_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
AM_M <- RunPCA(AM_M, features = VariableFeatures(object = AM_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
AM_M <- RunHarmony(AM_M, group.by.vars = "orig.ident")
AM_M <- RunUMAP(AM_M, reduction = "harmony", dims = 1:20)
AM_M <- RunTSNE(AM_M, reduction = "harmony", dims = 1:20)
names(AM_M@reductions)

###Cluster the cells.
AM_M <- FindNeighbors(AM_M, reduction = "harmony", dims = 1:20)
saveRDS(AM_M, file = "AM_Myeloid_before.rds")

#Determine resolution.
test <- AM_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_AM_Myeloid.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
AM_M <- FindClusters(AM_M, resolution = 0.1)
saveRDS(AM_M, file = "AM_Myeloid_res0.1.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="AM_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(AM_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="AM_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(AM_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="AM_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(AM_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="AM_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(AM_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="AM_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(AM_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="AM_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(AM_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="AM_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(AM_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="AM_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(AM_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="AM_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(AM_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="AM_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(AM_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="AM_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(AM_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="AM_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(AM_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="AM_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(AM_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="AM_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(AM_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="AM_Myeloid_features.pdf",width=10, height=14)
DotPlot(AM_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(4)
cDC2 = c(1)
C1QC_Macrophages = c(2)
LAMP3_DCs = c(5)
Langerhans_cells = c(8)
unknown_Myeloid = c(0, 3, 6, 7)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         C1QC_Macrophages,
                         LAMP3_DCs,
                         Langerhans_cells,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("Langerhans_cells",length(Langerhans_cells)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
)
AM_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AM_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(AM_M$celltype)
table(AM_M$seurat_clusters)

###Rename.
Idents(AM_M) <- AM_M$celltype
AM_M <- RenameIdents(AM_M, "LAMP3_DCs" = "LAMP3+ DCs", "C1QC_Macrophages" = "C1QC+ Macrophages", "Langerhans_cells" = "Langerhans cells", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(AM_M))
AM_M$celltype <- Idents(AM_M)
table(AM_M$celltype)
saveRDS(AM_M, file = "AM_Myeloid_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
AM_M <- subset(AM_M, subset = celltype %in% c("LAMP3+ DCs", "C1QC+ Macrophages", "Langerhans cells", "cDC2", "cDC1"))
AM_M$seurat_clusters <- as.factor(as.character(AM_M$seurat_clusters))
AM_M$celltype <- as.factor(as.character(AM_M$celltype))
AM_M$orig.ident <- as.factor(as.character(AM_M$orig.ident))

###Ranking.
Idents(AM_M) <- AM_M$celltype
AM_M <- RenameIdents(AM_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "Langerhans cells" = "Langerhans cells", "C1QC+ Macrophages" = "C1QC+ Macrophages")
table(Idents(AM_M))
AM_M$celltype <- Idents(AM_M)
table(AM_M$celltype)
saveRDS(AM_M, file = "AM_Myeloid_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6M: left.
pdf(file="AM_Myeloid_umap_celltype.pdf",width=6.3, height=4)
DimPlot(AM_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6M: right.
Idents(AM_M) <- 'celltype'
My_levels <- c("cDC1", "cDC2", "LAMP3+ DCs", "Langerhans cells", "C1QC+ Macrophages")
Idents(AM_M) <- factor(Idents(AM_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "C1QC", "CD207", "CD1A", "CCR7", "LAMP3", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="AM_Myeloid_M_features_res0.1.pdf",width=8, height=7.5)
DotPlot(AM_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6N.
###Cluster the T cells.
AM_T <- subset(AM, subset = Celltype %in% c("T cells"))
AM_T$seurat_clusters <- as.factor(as.character(AM_T$seurat_clusters))
AM_T$Celltype <- as.factor(as.character(AM_T$Celltype))
AM_T$orig.ident <- as.factor(as.character(AM_T$orig.ident))
AM_T$Sample <- as.factor(as.character(AM_T$Sample))
table(AM_T$Celltype)
saveRDS(AM_T, file = "AM_T_cells.rds")

###Normalizing the data.
AM_T <- NormalizeData(AM_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
AM_T <- FindVariableFeatures(AM_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(AM_T)
AM_T <- ScaleData(AM_T, features = all.genes)
AM_T <- ScaleData(AM_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
AM_T <- RunPCA(AM_T, features = VariableFeatures(object = AM_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
AM_T <- RunHarmony(AM_T, group.by.vars = "orig.ident")
AM_T <- RunUMAP(AM_T, reduction = "harmony", dims = 1:20)
AM_T <- RunTSNE(AM_T, reduction = "harmony", dims = 1:20)
names(AM_T@reductions)

###Cluster the cells.
AM_T <- FindNeighbors(AM_T, reduction = "harmony", dims = 1:20)
saveRDS(AM_T, file = "AM_T_cells_before.rds")

#Determine resolution.
test <- AM_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_AM_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
AM_T <- FindClusters(AM_T, resolution = 0.1)
saveRDS(AM_T, file = "AM_T_cells_res0.1.rds")

###classic celltype marker expression across different clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="AM_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(AM_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="AM_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(AM_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="AM_T_NKT.pdf",width=8, height=5)
DotPlot(AM_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="AM_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(AM_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4", "CD44", "CD69", "IL2RA")
pdf(file="AM_T_naive.pdf",width=8, height=4)
DotPlot(AM_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="AM_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(AM_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="AM_T_Exhaustion.pdf",width=8, height=4)
DotPlot(AM_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="AM_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(AM_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="AM_T_NK_marker.pdf",width=8, height=5)
DotPlot(AM_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "SELL", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="AM_T_cells_T_features.pdf",width=8, height=10)
DotPlot(AM_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(2)
CD4_T_cells = c(5, 7)
CD8_T_cells = c(1)
unknown_T_cells = c(0, 3, 4, 6)

current.cluster.ids <- c(Tregs,
                         CD4_T_cells,
                         CD8_T_cells,
                         unknown_T_cells)

new.cluster.ids <- c(rep("Tregs",length(Tregs)),
                     rep("CD4_T_cells",length(CD4_T_cells)),
                     rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)
AM_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AM_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(AM_T$seurat_clusters)
table(AM_T$celltype)

###Rename.
Idents(AM_T) <- AM_T$celltype
AM_T <- RenameIdents(AM_T, "CD8_T_cells" = "CD8+ T cells", "CD4_T_cells" = "CD4+ T cells", "unknown_T_cells" = "unknown T cells")
table(Idents(AM_T))
AM_T$celltype <- Idents(AM_T)
table(AM_T$celltype)
saveRDS(AM_T, file = "AM_T_cells_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
AM_T <- subset(AM_T, subset = celltype %in% c("Tregs", "CD4+ T cells", "CD8+ T cells"))
AM_T$seurat_clusters <- as.factor(as.character(AM_T$seurat_clusters))
AM_T$celltype <- as.factor(as.character(AM_T$celltype))
AM_T$orig.ident <- as.factor(as.character(AM_T$orig.ident))

###Ranking.
Idents(AM_T) <- AM_T$celltype
AM_T <- RenameIdents(AM_T, "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(AM_T))
AM_T$celltype <- Idents(AM_T)
table(AM_T$celltype)
saveRDS(AM_T, file = "AM_T_cells_res0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 6N: left.
pdf(file="AM_T_cells_umap_celltype.pdf",width=5.5, height=4)
DimPlot(AM_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6N: right.
Idents(AM_T) <- 'celltype'
My_levels <- c("Tregs", "CD4+ T cells", "CD8+ T cells")
Idents(AM_T) <- factor(Idents(AM_T), levels= My_levels)
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3")
pdf(file="AM_T_cells_T_features_res0.1_annotation.pdf",width=8, height=5)
DotPlot(AM_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
