###codes for Supplementary Figure 6O-P.

#CM.
#1 Supplementary Figure 6O.
###Cluster the Myeloid cells.
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM <- readRDS(file = "CM_resolution0.2_annotation2.rds") #The detailed code for CM_resolution0.2_annotation2.rds is at line 246 of the Supplementary Figure 5K-L.R file.
CM_M <- subset(CM, subset = Celltype %in% c("Myeloid cells"))
CM_M$seurat_clusters <- as.factor(as.character(CM_M$seurat_clusters))
CM_M$Celltype <- as.factor(as.character(CM_M$Celltype))
CM_M$orig.ident <- as.factor(as.character(CM_M$orig.ident))
CM_M$Sample <- as.factor(as.character(CM_M$Sample))
table(CM_M$Celltype)
saveRDS(CM_M, file = "CM_Myeloid.rds")

###Normalizing the data.
CM_M <- NormalizeData(CM_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CM_M <- FindVariableFeatures(CM_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CM_M)
CM_M <- ScaleData(CM_M, features = all.genes)
CM_M <- ScaleData(CM_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CM_M <- RunPCA(CM_M, features = VariableFeatures(object = CM_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CM_M <- RunHarmony(CM_M, group.by.vars = "orig.ident")
CM_M <- RunUMAP(CM_M, reduction = "harmony", dims = 1:20)
CM_M <- RunTSNE(CM_M, reduction = "harmony", dims = 1:20)
names(CM_M@reductions)

###Cluster the cells.
CM_M <- FindNeighbors(CM_M, reduction = "harmony", dims = 1:20)
saveRDS(CM_M, file = "CM_Myeloid_before.rds")

#Determine resolution.
test <- CM_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_CM_Myeloid.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.4.
CM_M <- FindClusters(CM_M, resolution = 0.4)
saveRDS(CM_M, file = "CM_Myeloid_res0.4.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="CM_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(CM_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="CM_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(CM_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="CM_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(CM_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="CM_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(CM_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="CM_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(CM_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="CM_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(CM_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="CM_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(CM_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="CM_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(CM_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="CM_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(CM_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="CM_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(CM_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="CM_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(CM_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="CM_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(CM_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="CM_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(CM_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="CM_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(CM_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="CM_Myeloid_features.pdf",width=10, height=14)
DotPlot(CM_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(5)
cDC2 = c(2)
C1QC_Macrophages = c(1)
pDCs = c(6)
Langerhans_cells = c(3)
unknown_Myeloid = c(0, 4)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         C1QC_Macrophages,
                         pDCs,
                         Langerhans_cells,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("pDCs",length(pDCs)),
                     rep("Langerhans_cells",length(Langerhans_cells)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
)

CM_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(CM_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CM_M$celltype)
table(CM_M$seurat_clusters)

###Rename.
Idents(CM_M) <- CM_M$celltype
CM_M <- RenameIdents(CM_M, "C1QC_Macrophages" = "C1QC+ Macrophages", "Langerhans_cells" = "Langerhans cells", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(CM_M))
CM_M$celltype <- Idents(CM_M)
table(CM_M$celltype)
saveRDS(CM_M, file = "CM_Myeloid_res0.4_annotation.rds")

###Exclude unknown cell subpopulations.
CM_M <- subset(CM_M, subset = celltype %in% c("pDCs", "C1QC+ Macrophages", "Langerhans cells", "cDC2", "cDC1"))
CM_M$seurat_clusters <- as.factor(as.character(CM_M$seurat_clusters))
CM_M$celltype <- as.factor(as.character(CM_M$celltype))
CM_M$orig.ident <- as.factor(as.character(CM_M$orig.ident))

###Ranking.
Idents(CM_M) <- CM_M$celltype
CM_M <- RenameIdents(CM_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "pDCs" = "pDCs", "Langerhans cells" = "Langerhans cells", "C1QC+ Macrophages" = "C1QC+ Macrophages")
table(Idents(CM_M))
CM_M$celltype <- Idents(CM_M)
table(CM_M$celltype)
saveRDS(CM_M, file = "CM_Myeloid_res0.4_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6O: left.
pdf(file="CM_Myeloid_umap_celltype.pdf",width=6.3, height=4)
DimPlot(CM_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6O: right.
Idents(CM_M) <- 'celltype'
My_levels <- c("cDC1", "cDC2", "pDCs", "Langerhans cells", "C1QC+ Macrophages")
Idents(CM_M) <- factor(Idents(CM_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "C1QC", "CD207", "CD1A", "PTGDS", "LILRA4", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="CM_Myeloid_M_features_res0.4.pdf",width=8, height=7.5)
DotPlot(CM_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6P.
###Cluster the T cells.
CM_T <- subset(CM, subset = Celltype %in% c("T cells"))
CM_T$seurat_clusters <- as.factor(as.character(CM_T$seurat_clusters))
CM_T$Celltype <- as.factor(as.character(CM_T$Celltype))
CM_T$orig.ident <- as.factor(as.character(CM_T$orig.ident))
CM_T$Sample <- as.factor(as.character(CM_T$Sample))
table(CM_T$Celltype)
saveRDS(CM_T, file = "CM_T_cells.rds")

###Normalizing the data.
CM_T <- NormalizeData(CM_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CM_T <- FindVariableFeatures(CM_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CM_T)
CM_T <- ScaleData(CM_T, features = all.genes)
CM_T <- ScaleData(CM_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CM_T <- RunPCA(CM_T, features = VariableFeatures(object = CM_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CM_T <- RunHarmony(CM_T, group.by.vars = "orig.ident")
CM_T <- RunUMAP(CM_T, reduction = "harmony", dims = 1:20)
CM_T <- RunTSNE(CM_T, reduction = "harmony", dims = 1:20)
names(CM_T@reductions)

###Cluster the cells.
CM_T <- FindNeighbors(CM_T, reduction = "harmony", dims = 1:20)
saveRDS(CM_T, file = "CM_T_cells_before.rds")

#Determine resolution.
test <- CM_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_CM_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
CM_T <- FindClusters(CM_T, resolution = 0.2)
saveRDS(CM_T, file = "CM_T_cells_res0.2.rds")

###classic celltype marker expression across different clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="CM_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(CM_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="CM_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(CM_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="CM_T_NKT.pdf",width=8, height=5)
DotPlot(CM_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="CM_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(CM_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4", "CD44", "CD69", "IL2RA")
pdf(file="CM_T_naive.pdf",width=8, height=4)
DotPlot(CM_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="CM_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(CM_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="CM_T_Exhaustion.pdf",width=8, height=4)
DotPlot(CM_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="CM_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(CM_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="CM_T_NK_marker.pdf",width=8, height=5)
DotPlot(CM_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "SELL", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="CM_T_cells_T_features.pdf",width=8, height=10)
DotPlot(CM_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(2, 3)
NK_cells = c(0, 5)
CD8_T_cells = c(6)
unknown_T_cells = c(1, 4, 7, 8)
current.cluster.ids <- c(Tregs,
                         NK_cells,
                         CD8_T_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("Tregs",length(Tregs)),
                     rep("NK_cells",length(NK_cells)),
                     rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

CM_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(CM_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CM_T$seurat_clusters)
table(CM_T$celltype)

###Rename.
Idents(CM_T) <- CM_T$celltype
CM_T <- RenameIdents(CM_T, "CD8_T_cells" = "CD8+ T cells", "NK_cells" = "NK cells", "unknown_T_cells" = "unknown T cells")
table(Idents(CM_T))
CM_T$celltype <- Idents(CM_T)
table(CM_T$celltype)
saveRDS(CM_T, file = "CM_T_cells_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
CM_T <- subset(CM_T, subset = celltype %in% c("CD8+ T cells", "Tregs", "NK cells"))
CM_T$seurat_clusters <- as.factor(as.character(CM_T$seurat_clusters))
CM_T$celltype <- as.factor(as.character(CM_T$celltype))

###Ranking.
Idents(CM_T) <- CM_T$celltype
CM_T <- RenameIdents(CM_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD8+ T cells" = "CD8+ T cells")
table(Idents(CM_T))
CM_T$celltype <- Idents(CM_T)
table(CM_T$celltype)
saveRDS(CM_T, file = "CM_T_cells_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6P: left.
pdf(file="CM_T_cells_umap_celltype.pdf",width=5.6, height=4)
DimPlot(CM_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6P: right.
Idents(CM_T) <- 'celltype'
My_levels <- c("NK cells", "Tregs", "CD8+ T cells")
Idents(CM_T) <- factor(Idents(CM_T), levels= My_levels)
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="CM_T_cells_T_features_res0.2_annotation.pdf",width=8, height=5)
DotPlot(CM_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
