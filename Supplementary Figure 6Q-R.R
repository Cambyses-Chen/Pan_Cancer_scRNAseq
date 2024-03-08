###codes for Supplementary Figure 6Q-R.

#RCC.
#1 Supplementary Figure 6Q.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
library(Seurat)
library(scuttle)
library(ggplot2)
library(scales)
library(ggsci)
RCC <- readRDS(file = "RCC_resolution0.1_annotation2.rds") #The detailed code for RCC_resolution0.1_annotation2.rds is at line 288 of the Supplementary Figure 5M-N.R file.
RCC_M <- subset(RCC, subset = Celltype %in% c("Myeloid cells"))
RCC_M$seurat_clusters <- as.factor(as.character(RCC_M$seurat_clusters))
RCC_M$Celltype <- as.factor(as.character(RCC_M$Celltype))
RCC_M$orig.ident <- as.factor(as.character(RCC_M$orig.ident))
saveRDS(RCC_M, file = "RCC_Myeloid.rds")

###Normalizing the data.
RCC_M <- NormalizeData(RCC_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
RCC_M <- FindVariableFeatures(RCC_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(RCC_M)
RCC_M <- ScaleData(RCC_M, features = all.genes)
RCC_M <- ScaleData(RCC_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
RCC_M <- RunPCA(RCC_M, features = VariableFeatures(object = RCC_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
RCC_M <- RunHarmony(RCC_M, group.by.vars = "orig.ident")
RCC_M <- RunUMAP(RCC_M, reduction = "harmony", dims = 1:20)
RCC_M <- RunTSNE(RCC_M, reduction = "harmony", dims = 1:20)
names(RCC_M@reductions)

###Cluster the cells.
RCC_M <- FindNeighbors(RCC_M, reduction = "harmony", dims = 1:20)
saveRDS(RCC_M, file = "RCC_Myeloid_before.rds")

#Determine resolution.
test <- RCC_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_RCC_M.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
RCC_M <- FindClusters(RCC_M, resolution = 0.1)
saveRDS(RCC_M, file = "RCC_Myeloid_res0.1.rds")

###classic celltype marker in clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="RCC_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(RCC_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="RCC_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(RCC_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="RCC_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(RCC_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="RCC_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(RCC_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="RCC_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(RCC_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="RCC_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(RCC_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="RCC_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(RCC_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="RCC_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(RCC_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="RCC_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(RCC_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="RCC_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(RCC_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="RCC_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(RCC_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="RCC_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(RCC_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="RCC_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(RCC_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="RCC_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(RCC_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="RCC_Myeloid_features.pdf",width=10, height=14)
DotPlot(RCC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(2)
cDC2 = c(3)
pDCs = c(4)
SPP1_Macrophages = c(5)
VCAN_Macrophages = c(1)
C1QC_Macrophages = c(0)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         pDCs,
                         SPP1_Macrophages,
                         VCAN_Macrophages,
                         C1QC_Macrophages)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("pDCs",length(pDCs)),
                     rep("SPP1_Macrophages",length(SPP1_Macrophages)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages))
                     
)

RCC_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(RCC_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(RCC_M$celltype)

###Rename.
Idents(RCC_M) <- RCC_M$celltype
RCC_M <- RenameIdents(RCC_M, "SPP1_Macrophages" = "SPP1+ Macrophages", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages")
table(Idents(RCC_M))
RCC_M$celltype <- Idents(RCC_M)
table(RCC_M$celltype)
saveRDS(RCC_M, file = "RCC_Myeloid_res0.1_annotation.rds")

###Ranking.
Idents(RCC_M) <- RCC_M$celltype
RCC_M <- RenameIdents(RCC_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "pDCs" = "pDCs", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(RCC_M))
RCC_M$celltype <- Idents(RCC_M)
table(RCC_M$celltype)
saveRDS(RCC_M, file = "RCC_Myeloid_res0.1_annotation2.rds")

###Visualization.
#UMAP plot for Supplementary Figure 6Q: left.
pdf(file="RCC_Myeloid_umap_celltype.pdf",width=6, height=4)
DimPlot(RCC_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6Q: right.
Idents(RCC_M) <- 'celltype'
My_levels <- c( "cDC1", "cDC2", "pDCs", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages")
Idents(RCC_M) <- factor(Idents(RCC_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "SPP1", "C1QC", "PTGDS", "LILRA4", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="RCC_Myeloid_Myeloid_features.pdf",width=8, height=7.5)
DotPlot(RCC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6R.
###Cluster the T cells.
RCC_T <- subset(RCC, subset = Celltype %in% c("T cells"))
RCC_T$seurat_clusters <- as.factor(as.character(RCC_T$seurat_clusters))
RCC_T$Celltype <- as.factor(as.character(RCC_T$Celltype))
RCC_T$orig.ident <- as.factor(as.character(RCC_T$orig.ident))
table(RCC_T$Celltype)
saveRDS(RCC_T, file = "RCC_T_cells.rds")

###Normalizing the data.
RCC_T <- NormalizeData(RCC_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
RCC_T <- FindVariableFeatures(RCC_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(RCC_T)
RCC_T <- ScaleData(RCC_T, features = all.genes)
RCC_T <- ScaleData(RCC_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
RCC_T <- RunPCA(RCC_T, features = VariableFeatures(object = RCC_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
RCC_T <- RunHarmony(RCC_T, group.by.vars = "orig.ident")
RCC_T <- RunUMAP(RCC_T, reduction = "harmony", dims = 1:20)
RCC_T <- RunTSNE(RCC_T, reduction = "harmony", dims = 1:20)
names(RCC_T@reductions)

###Cluster the cells.
RCC_T <- FindNeighbors(RCC_T, reduction = "harmony", dims = 1:20)
saveRDS(RCC_T, file = "RCC_T_cells_before.rds")

#Determine resolution.
test <- RCC_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_RCC_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
RCC_T <- FindClusters(RCC_T, resolution = 0.1)
saveRDS(RCC_T, file = "RCC_T_cells_res0.1.rds")

###classic celltype marker in clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="RCC_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(RCC_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="RCC_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(RCC_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="RCC_T_NKT.pdf",width=8, height=5)
DotPlot(RCC_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="RCC_T_CD4_Treg_0213.pdf",width=8, height=4)
DotPlot(RCC_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="RCC_T_naive.pdf",width=8, height=4)
DotPlot(RCC_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="RCC_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(RCC_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="RCC_T_Exhaustion.pdf",width=8, height=4)
DotPlot(RCC_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="RCC_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(RCC_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="RCC_T_NK_marker.pdf",width=8, height=5)
DotPlot(RCC_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="RCC_T_cells_T_features_res0.1_annotation.pdf",width=8, height=10)
DotPlot(RCC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
CD8_T_cells = c(0)
Tregs = c(4)
CD4_T_cells = c(1)
unknown_T_cells = c(2, 3)
current.cluster.ids <- c(CD8_T_cells,
                         Tregs,
                         CD4_T_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("Tregs",length(Tregs)),
                     rep("CD4_T_cells",length(CD4_T_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

RCC_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(RCC_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(RCC_T$celltype)
table(RCC_T$seurat_clusters)

###Rename.
Idents(RCC_T) <- RCC_T$celltype
RCC_T <- RenameIdents(RCC_T, "CD8_T_cells" = "CD8+ T cells", "CD4_T_cells" = "CD4+ T cells", "unknown_T_cells" = "unknown T cells")
table(Idents(RCC_T))
RCC_T$celltype <- Idents(RCC_T)
table(RCC_T$celltype)
saveRDS(RCC_T, file = "RCC_T_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
RCC_T <- subset(RCC_T, subset = celltype %in% c("CD8+ T cells", "CD4+ T cells", "Tregs"))
RCC_T$seurat_clusters <- as.factor(as.character(RCC_T$seurat_clusters))
RCC_T$celltype <- as.factor(as.character(RCC_T$celltype))
RCC_T$orig.ident <- as.factor(as.character(RCC_T$orig.ident))

###Ranking.
Idents(RCC_T) <- RCC_T$celltype
RCC_T <- RenameIdents(RCC_T, "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(RCC_T))
RCC_T$celltype <- Idents(RCC_T)
table(RCC_T$celltype)
saveRDS(RCC_T, file = "RCC_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6R: left.
pdf(file="RCC_T_cells_umap_celltype.pdf",width=5.7, height=4)
DimPlot(RCC_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6R: right.
Idents(RCC_T) <- RCC_T$celltype
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3")
pdf(file="RCC_T_cells_T_features_res0.1_annotation.pdf",width=8, height=5)
DotPlot(RCC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
