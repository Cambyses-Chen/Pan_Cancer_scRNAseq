###codes for Supplementary Figure 6C-D.

#OV.
#1 Supplementary Figure 6C.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
OV <- readRDS(file = "OV_resolution0.2_annotation2.rds") #The detailed code for OV_resolution0.2_annotation2.rds is at line 261 of the Figure 1G-J.R file.
OV_M <- subset(OV, subset = Celltype %in% c("Myeloid cells"))
OV_M$seurat_clusters <- as.factor(as.character(OV_M$seurat_clusters))
OV_M$Celltype <- as.factor(as.character(OV_M$Celltype))
OV_M$orig.ident <- as.factor(as.character(OV_M$orig.ident))
OV_M$Treatment <- as.factor(as.character(OV_M$Treatment))
table(OV_M$Celltype)
saveRDS(OV_M, file = "OV_Myeloid.rds")

###Normalizing the data.
OV_M <- NormalizeData(OV_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
OV_M <- FindVariableFeatures(OV_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(OV_M)
OV_M <- ScaleData(OV_M, features = all.genes)
OV_M <- ScaleData(OV_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
OV_M <- RunPCA(OV_M, features = VariableFeatures(object = OV_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
OV_M <- RunHarmony(OV_M, group.by.vars = "orig.ident")
OV_M <- RunUMAP(OV_M, reduction = "harmony", dims = 1:20)
OV_M <- RunTSNE(OV_M, reduction = "harmony", dims = 1:20)
names(OV_M@reductions)

###Cluster the cells.
OV_M <- FindNeighbors(OV_M, reduction = "harmony", dims = 1:20)
saveRDS(OV_M, file = "OV_Myeloid_before.rds")

#Determine resolution.
test <- OV_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_OV_M.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.4.
OV_M <- FindClusters(OV_M, resolution = 0.4)
saveRDS(OV_M, file = "OV_Myeloid_res0.4.rds")

###classic celltype marker in clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="OV_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(OV_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="OV_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(OV_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="OV_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(OV_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="OV_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(OV_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="OV_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(OV_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="OV_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(OV_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="OV_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(OV_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="OV_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(OV_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="OV_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(OV_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="OV_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(OV_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="OV_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(OV_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="OV_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(OV_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="OV_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(OV_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="OV_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(OV_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="OV_Myeloid_features.pdf",width=10, height=14)
DotPlot(OV_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(11)
cDC2 = c(4)
pDCs = c(7)
SPP1_Macrophages = c(5, 8)
VCAN_Macrophages = c(1)
C1QC_Macrophages = c(2, 3, 6, 10)
unknown_Myeloid = c(0, 9, 12)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         pDCs,
                         SPP1_Macrophages,
                         VCAN_Macrophages,
                         C1QC_Macrophages,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("pDCs",length(pDCs)),
                     rep("SPP1_Macrophages",length(SPP1_Macrophages)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
                     
)

OV_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(OV_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(OV_M$celltype)
table(OV_M$seurat_clusters)

###Rename.
Idents(OV_M) <- OV_M$celltype
OV_M <- RenameIdents(OV_M, "SPP1_Macrophages" = "SPP1+ Macrophages", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(OV_M))
OV_M$celltype <- Idents(OV_M)
table(OV_M$celltype)
saveRDS(OV_M, file = "OV_Myeloid_res0.4_annotation.rds")

###Exclude unknown cell subpopulations.
OV_M <- subset(OV_M, subset = celltype %in% c("pDCs", "C1QC+ Macrophages", "VCAN+ Macrophages", "cDC1", "cDC2", "SPP1+ Macrophages"))
OV_M$seurat_clusters <- as.factor(as.character(OV_M$seurat_clusters))
OV_M$celltype <- as.factor(as.character(OV_M$celltype))
OV_M$orig.ident <- as.factor(as.character(OV_M$orig.ident))

###Ranking.
Idents(OV_M) <- OV_M$celltype
OV_M <- RenameIdents(OV_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "pDCs" = "pDCs", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(OV_M))
OV_M$celltype <- Idents(OV_M)
table(OV_M$celltype)
saveRDS(OV_M, file = "OV_Myeloid_res0.4_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6C: left.
pdf(file="OV_Myeloid_umap_celltype.pdf",width=6, height=4)
DimPlot(OV_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6C: right.
Idents(OV_M) <- 'celltype'
My_levels <- c( "cDC1", "cDC2", "pDCs", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages")
Idents(OV_M) <- factor(Idents(OV_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "SPP1", "C1QC", "PTGDS", "LILRA4", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="OV_Myeloid_Myeloid_features.pdf",width=8, height=8)
DotPlot(OV_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6D.
###Cluster the T cells.
OV_T <- subset(OV, subset = Celltype %in% c("T cells"))
OV_T$seurat_clusters <- as.factor(as.character(OV_T$seurat_clusters))
OV_T$Celltype <- as.factor(as.character(OV_T$Celltype))
OV_T$orig.ident <- as.factor(as.character(OV_T$orig.ident))
OV_T$Treatment <- as.factor(as.character(OV_T$Treatment))
saveRDS(OV_T, file = "OV_T.rds")

###Normalizing the data.
OV_T <- NormalizeData(OV_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
OV_T <- FindVariableFeatures(OV_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(OV_T)
OV_T <- ScaleData(OV_T, features = all.genes)
OV_T <- ScaleData(OV_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
OV_T <- RunPCA(OV_T, features = VariableFeatures(object = OV_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
OV_T <- RunHarmony(OV_T, group.by.vars = "orig.ident")
OV_T <- RunUMAP(OV_T, reduction = "harmony", dims = 1:20)
OV_T <- RunTSNE(OV_T, reduction = "harmony", dims = 1:20)
names(OV_T@reductions)

###Cluster the cells.
OV_T <- FindNeighbors(OV_T, reduction = "harmony", dims = 1:20)
saveRDS(OV_T, file = "OV_T_cells_before.rds")

#Determine resolution.
test <- OV_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_OV_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
OV_T <- FindClusters(OV_T, resolution = 0.1)
saveRDS(OV_T, file = "OV_T_cells_res0.1.rds")

###classic celltype marker in clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="OV_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(OV_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="OV_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(OV_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="OV_T_NKT.pdf",width=8, height=5)
DotPlot(OV_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="OV_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(OV_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="OV_T_naive.pdf",width=8, height=4)
DotPlot(OV_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="OV_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(OV_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="OV_T_Exhaustion.pdf",width=8, height=4)
DotPlot(OV_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="OV_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(OV_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="OV_T_NK_marker.pdf",width=8, height=5)
DotPlot(OV_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "SELL", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="OV_T_cells_T_features.pdf",width=8, height=10)
DotPlot(OV_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(1)
CD8_T_cells = c(0)
NK_cells = c(2)
unknown_T_cells = c(3, 4, 5, 6)
current.cluster.ids <- c(Tregs,
                         CD8_T_cells,
                         NK_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("Tregs",length(Tregs)),
                     rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

OV_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(OV_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(OV_T$seurat_clusters)
table(OV_T$celltype)

###Rename.
Idents(OV_T) <- OV_T$celltype
OV_T <- RenameIdents(OV_T, "CD8_T_cells" = "CD8+ T cells", "NK_cells" = "NK cells", "unknown_T_cells" = "unknown T cells")
table(Idents(OV_T))
OV_T$celltype <- Idents(OV_T)
table(OV_T$celltype)
saveRDS(OV_T, file = "OV_T_cells_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
OV_T <- subset(OV_T, subset = celltype %in% c("CD8+ T cells", "NK cells", "Tregs"))
OV_T$seurat_clusters <- as.factor(as.character(OV_T$seurat_clusters))
OV_T$celltype <- as.factor(as.character(OV_T$celltype))

###Ranking.
Idents(OV_T) <- OV_T$celltype
OV_T <- RenameIdents(OV_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD8+ T cells" = "CD8+ T cells")
table(Idents(OV_T))
OV_T$celltype <- Idents(OV_T)
table(OV_T$celltype)
saveRDS(OV_T, file = "OV_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6D: left.
pdf(file="OV_T_cells_umap_celltype.pdf",width=5.5, height=4)
DimPlot(OV_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6D: right.
Idents(OV_T) <- 'celltype'
My_levels <- c("NK cells", "Tregs", "CD8+ T cells")
Idents(OV_T) <- factor(Idents(OV_T), levels= My_levels)
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="OV_T_cells_T_features_res0.1_annotation.pdf",width=8, height=4.5)
DotPlot(OV_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
