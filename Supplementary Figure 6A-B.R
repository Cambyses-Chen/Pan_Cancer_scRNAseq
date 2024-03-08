###codes for Supplementary Figure 6A-B.

#HNSCC.
#1 Supplementary Figure 6A.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
HNSCC <- readRDS(file = "HNSCC_resolution0.3_annotation_sub.rds") #The detailed code for HNSCC_resolution0.3_annotation_sub.rds is at line 273 of the Figure 1C-F.R file.
HNSCC_M <- subset(HNSCC, subset = Celltype %in% c("Myeloid cells"))
HNSCC_M$seurat_clusters <- as.factor(as.character(HNSCC_M$seurat_clusters))
HNSCC_M$Celltype <- as.factor(as.character(HNSCC_M$Celltype))
HNSCC_M$orig.ident <- as.factor(as.character(HNSCC_M$orig.ident))
saveRDS(HNSCC_M, file = "HNSCC_Myeloid.rds")

###Normalizing the data.
HNSCC_M <- NormalizeData(HNSCC_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
HNSCC_M <- FindVariableFeatures(HNSCC_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(HNSCC_M)
HNSCC_M <- ScaleData(HNSCC_M, features = all.genes)
HNSCC_M <- ScaleData(HNSCC_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
HNSCC_M <- RunPCA(HNSCC_M, features = VariableFeatures(object = HNSCC_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
HNSCC_M <- RunHarmony(HNSCC_M, group.by.vars = "orig.ident")
HNSCC_M <- RunUMAP(HNSCC_M, reduction = "harmony", dims = 1:20)
HNSCC_M <- RunTSNE(HNSCC_M, reduction = "harmony", dims = 1:20)
names(HNSCC_M@reductions)

###Cluster the cells.
HNSCC_M <- FindNeighbors(HNSCC_M, reduction = "harmony", dims = 1:20)
saveRDS(HNSCC_M, file = "HNSCC_Myeloid_before.rds")

#Determine resolution.
test <- HNSCC_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_HNSCC_M.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.4.
HNSCC_M <- FindClusters(HNSCC_M, resolution = 0.4)
saveRDS(HNSCC_M, file = "HNSCC_Myeloid_res0.4.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="HNSCC_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(HNSCC_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="HNSCC_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(HNSCC_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="HNSCC_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(HNSCC_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="HNSCC_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(HNSCC_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="HNSCC_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(HNSCC_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="HNSCC_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(HNSCC_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="HNSCC_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(HNSCC_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="HNSCC_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(HNSCC_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="HNSCC_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(HNSCC_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="HNSCC_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(HNSCC_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="HNSCC_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(HNSCC_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="HNSCC_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(HNSCC_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="HNSCC_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(HNSCC_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="HNSCC_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(HNSCC_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="HNSCC_Myeloid_features.pdf",width=10, height=14)
DotPlot(HNSCC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(6)
cDC2 = c(2, 8)
LAMP3_DCs = c(5)
SPP1_Macrophages = c(0, 4)
LYVE1_Macrophages = c(3)
C1QC_Macrophages = c(7, 10)
INHBA_Macrophages = c(9)
unknown_Myeloid = c(1, 11, 12, 13)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         LAMP3_DCs,
                         SPP1_Macrophages,
                         LYVE1_Macrophages,
                         C1QC_Macrophages,
                         INHBA_Macrophages,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("SPP1_Macrophages",length(SPP1_Macrophages)),
                     rep("LYVE1_Macrophages",length(LYVE1_Macrophages)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("INHBA_Macrophages",length(INHBA_Macrophages)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
                     
)

HNSCC_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(HNSCC_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(HNSCC_M$celltype)
table(HNSCC_M$seurat_clusters)

###Rename.
Idents(HNSCC_M) <- HNSCC_M$celltype
HNSCC_M <- RenameIdents(HNSCC_M, "SPP1_Macrophages" = "SPP1+ Macrophages", "C1QC_Macrophages" = "C1QC+ Macrophages", "INHBA_Macrophages" = "INHBA+ Macrophages", "LYVE1_Macrophages" = "LYVE1+ Macrophages", "LAMP3_DCs" = "LAMP3+ DCs", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(HNSCC_M))
HNSCC_M$celltype <- Idents(HNSCC_M)
table(HNSCC_M$celltype)
saveRDS(HNSCC_M, file = "HNSCC_Myeloid_res0.4_annotation.rds")

###Exclude unknown cell subpopulations.
HNSCC_M <- subset(HNSCC_M, subset = celltype %in% c("SPP1+ Macrophages", "C1QC+ Macrophages", "INHBA+ Macrophages", "LYVE1+ Macrophages", "LAMP3+ DCs", "cDC2", "cDC1"))
HNSCC_M$seurat_clusters <- as.factor(as.character(HNSCC_M$seurat_clusters))
HNSCC_M$celltype <- as.factor(as.character(HNSCC_M$celltype))
HNSCC_M$orig.ident <- as.factor(as.character(HNSCC_M$orig.ident))
HNSCC_M$Sample <- as.factor(as.character(HNSCC_M$Sample))

###Ranking.
Idents(HNSCC_M) <- HNSCC_M$celltype
HNSCC_M <- RenameIdents(HNSCC_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "LYVE1+ Macrophages" = "LYVE1+ Macrophages", "INHBA+ Macrophages" = "INHBA+ Macrophages")
table(Idents(HNSCC_M))
HNSCC_M$celltype <- Idents(HNSCC_M)
table(HNSCC_M$celltype)
saveRDS(HNSCC_M, file = "HNSCC_Myeloid_res0.4_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6A: left.
pdf(file="HNSCC_Myeloid_umap_celltype.pdf",width=6.2, height=4)
DimPlot(HNSCC_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6A: right.
Idents(HNSCC_M) <- 'celltype'
My_levels <- c( "cDC1", "cDC2", "LAMP3+ DCs", "C1QC+ Macrophages", "SPP1+ Macrophages", "LYVE1+ Macrophages", "INHBA+ Macrophages")
Idents(HNSCC_M) <- factor(Idents(HNSCC_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "INHBA", "LYVE1", "SPP1", "C1QC", "CCR7", "LAMP3", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="HNSCC_Myeloid_Myeloid_features.pdf",width=8, height=8)
DotPlot(HNSCC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6B.
###Cluster the T cells.
HNSCC_T <- subset(HNSCC, subset = Celltype %in% c("T cells"))
HNSCC_T$seurat_clusters <- as.factor(as.character(HNSCC_T$seurat_clusters))
HNSCC_T$Celltype <- as.factor(as.character(HNSCC_T$Celltype))
HNSCC_T$orig.ident <- as.factor(as.character(HNSCC_T$orig.ident))
table(HNSCC_T$Celltype)
saveRDS(HNSCC_T, file = "HNSCC_T_cells.rds")

###Normalizing the data.
HNSCC_T <- NormalizeData(HNSCC_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
HNSCC_T <- FindVariableFeatures(HNSCC_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(HNSCC_T)
HNSCC_T <- ScaleData(HNSCC_T, features = all.genes)
HNSCC_T <- ScaleData(HNSCC_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
HNSCC_T <- RunPCA(HNSCC_T, features = VariableFeatures(object = HNSCC_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
HNSCC_T <- RunHarmony(HNSCC_T, group.by.vars = "orig.ident")
HNSCC_T <- RunUMAP(HNSCC_T, reduction = "harmony", dims = 1:20)
HNSCC_T <- RunTSNE(HNSCC_T, reduction = "harmony", dims = 1:20)
names(HNSCC_T@reductions)

###Cluster the cells.
HNSCC_T <- FindNeighbors(HNSCC_T, reduction = "harmony", dims = 1:20)
saveRDS(HNSCC_T, file = "HNSCC_T_cells_before.rds")

#Determine resolution.
test <- HNSCC_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_HNSCC_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
HNSCC_T <- FindClusters(HNSCC_T, resolution = 0.1)
saveRDS(HNSCC_T, file = "HNSCC_T_cells_res0.1.rds")

###classic celltype marker in clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="HNSCC_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="HNSCC_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="HNSCC_T_NKT.pdf",width=8, height=5)
DotPlot(HNSCC_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="HNSCC_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="HNSCC_T_naive.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="HNSCC_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(HNSCC_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="HNSCC_T_Exhaustion.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="HNSCC_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(HNSCC_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="HNSCC_T_NK_marker.pdf",width=8, height=5)
DotPlot(HNSCC_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="HNSCC_T_cells_T_features_res0.2_annotation.pdf",width=8, height=10)
DotPlot(HNSCC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
CD8_T_cells = c(0)
Tregs = c(2)
NK_cells = c(4)
unknown_T_cells = c(1, 3, 5, 6)
current.cluster.ids <- c(CD8_T_cells,
                         Tregs,
                         NK_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("Tregs",length(Tregs)),
                     rep("NK_cells",length(NK_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

HNSCC_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(HNSCC_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(HNSCC_T$celltype)
table(HNSCC_T$seurat_clusters)

###Rename.
Idents(HNSCC_T) <- HNSCC_T$celltype
HNSCC_T <- RenameIdents(HNSCC_T, "CD8_T_cells" = "CD8+ T cells", "NK_cells" = "NK cells", "unknown_T_cells" = "unknown T cells")
table(Idents(HNSCC_T))
HNSCC_T$celltype <- Idents(HNSCC_T)
table(HNSCC_T$celltype)
saveRDS(HNSCC_T, file = "HNSCC_T_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
HNSCC_T <- subset(HNSCC_T, subset = celltype %in% c("CD8+ T cells", "NK cells", "Tregs"))
HNSCC_T$seurat_clusters <- as.factor(as.character(HNSCC_T$seurat_clusters))
HNSCC_T$celltype <- as.factor(as.character(HNSCC_T$celltype))
HNSCC_T$orig.ident <- as.factor(as.character(HNSCC_T$orig.ident))

###Ranking.
Idents(HNSCC_T) <- HNSCC_T$celltype
HNSCC_T <- RenameIdents(HNSCC_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD8+ T cells" = "CD8+ T cells")
table(Idents(HNSCC_T))
HNSCC_T$celltype <- Idents(HNSCC_T)
table(HNSCC_T$celltype)
saveRDS(HNSCC_T, file = "HNSCC_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6B: left.
pdf(file="HNSCC_T_cells_umap_celltype.pdf",width=5.6, height=4)
DimPlot(HNSCC_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6B: right.
Idents(HNSCC_T) <- HNSCC_T$celltype
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="HNSCC_T_cells_T_features_res0.1_annotation.pdf",width=8, height=5)
DotPlot(HNSCC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
