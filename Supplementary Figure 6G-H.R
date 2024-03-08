###codes for Supplementary Figure 6G-H.

#CC.
#1 Supplementary Figure 6G.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
Cervical <- readRDS(file = "Cervical_resolution0.3_annotation_2.rds") #The detailed code for Cervical_resolution0.3_annotation_2.rds is at line 268 of the Supplementary Figure 5C-D.R file.
Cervical_M <- subset(Cervical, subset = Celltype %in% c("Myeloid cells"))
Cervical_M$seurat_clusters <- as.factor(as.character(Cervical_M$seurat_clusters))
Cervical_M$Celltype <- as.factor(as.character(Cervical_M$Celltype))
Cervical_M$orig.ident <- as.factor(as.character(Cervical_M$orig.ident))
saveRDS(Cervical_M, file = "Cervical_Myeloid.rds")

###Normalizing the data.
Cervical_M <- NormalizeData(Cervical_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
Cervical_M <- FindVariableFeatures(Cervical_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(Cervical_M)
Cervical_M <- ScaleData(Cervical_M, features = all.genes)
Cervical_M <- ScaleData(Cervical_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
Cervical_M <- RunPCA(Cervical_M, features = VariableFeatures(object = Cervical_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
Cervical_M <- RunHarmony(Cervical_M, group.by.vars = "orig.ident")
Cervical_M <- RunUMAP(Cervical_M, reduction = "harmony", dims = 1:20)
Cervical_M <- RunTSNE(Cervical_M, reduction = "harmony", dims = 1:20)
names(Cervical_M@reductions)

###Cluster the cells.
Cervical_M <- FindNeighbors(Cervical_M, reduction = "harmony", dims = 1:20)
saveRDS(Cervical_M, file = "Cervical_Myeloid_before.rds")

#Determine resolution.
test <- Cervical_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_Cervical_M.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
Cervical_M <- FindClusters(Cervical_M, resolution = 0.2)
saveRDS(Cervical_M, file = "Cervical_Myeloid_res0.2.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="Cervical_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(Cervical_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="Cervical_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(Cervical_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="Cervical_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(Cervical_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="Cervical_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(Cervical_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="Cervical_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(Cervical_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="Cervical_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(Cervical_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="Cervical_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(Cervical_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="Cervical_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(Cervical_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="Cervical_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(Cervical_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="Cervical_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(Cervical_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="Cervical_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(Cervical_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="Cervical_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(Cervical_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="Cervical_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(Cervical_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="Cervical_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(Cervical_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="Cervical_Myeloid_features.pdf",width=10, height=14)
DotPlot(Cervical_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(7)
LAMP3_DCs = c(9)
C1QC_Macrophages = c(1)
SPP1_Macrophages = c(6)
VCAN_Macrophages = c(4)
Neutrophils = c(0, 2, 5, 10)
unknown_Myeloid = c(3, 8)
current.cluster.ids <- c(cDC1,
                         LAMP3_DCs,
                         C1QC_Macrophages,
                         SPP1_Macrophages,
                         VCAN_Macrophages,
                         Neutrophils,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("SPP1_Macrophages",length(SPP1_Macrophages)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("Neutrophils",length(Neutrophils)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
                     
)

Cervical_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Cervical_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Cervical_M$celltype)
table(Cervical_M$seurat_clusters)

###Rename.
Idents(Cervical_M) <- Cervical_M$celltype
Cervical_M <- RenameIdents(Cervical_M, "SPP1_Macrophages" = "SPP1+ Macrophages", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages", "LAMP3_DCs" = "LAMP3+ DCs", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(Cervical_M))
Cervical_M$celltype <- Idents(Cervical_M)
table(Cervical_M$celltype)
saveRDS(Cervical_M, file = "Cervical_Myeloid_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
Cervical_M <- subset(Cervical_M, subset = celltype %in% c("SPP1+ Macrophages", "C1QC+ Macrophages", "VCAN+ Macrophages", "Neutrophils", "LAMP3+ DCs", "cDC1"))
Cervical_M$seurat_clusters <- as.factor(as.character(Cervical_M$seurat_clusters))
Cervical_M$celltype <- as.factor(as.character(Cervical_M$celltype))
Cervical_M$orig.ident <- as.factor(as.character(Cervical_M$orig.ident))
Cervical_M$Sample <- as.factor(as.character(Cervical_M$Sample))

###Ranking.
Idents(Cervical_M) <- Cervical_M$celltype
Cervical_M <- RenameIdents(Cervical_M, "cDC1" = "cDC1", "LAMP3+ DCs" = "LAMP3+ DCs", "Neutrophils" = "Neutrophils", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(Cervical_M))
Cervical_M$celltype <- Idents(Cervical_M)
table(Cervical_M$celltype)
saveRDS(Cervical_M, file = "Cervical_Myeloid_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6G: left.
pdf(file="Cervical_Myeloid_umap_celltype.pdf",width=6.2, height=4)
DimPlot(Cervical_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6G: right.
Idents(Cervical_M) <- 'celltype'
My_levels <- c( "cDC1", "LAMP3+ DCs", "Neutrophils", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages")
Idents(Cervical_M) <- factor(Idents(Cervical_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "SPP1", "C1QC", "KATNBL1", "FCGR3B", "CCR7", "LAMP3", "CLEC9A", "CADM1")
pdf(file="Cervical_Myeloid_Myeloid_features.pdf",width=8, height=8)
DotPlot(Cervical_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6H.
###Cluster the T cells.
Cervical_T <- subset(Cervical, subset = Celltype %in% c("T cells"))
Cervical_T$seurat_clusters <- as.factor(as.character(Cervical_T$seurat_clusters))
Cervical_T$Celltype <- as.factor(as.character(Cervical_T$Celltype))
Cervical_T$orig.ident <- as.factor(as.character(Cervical_T$orig.ident))
table(Cervical_T$Celltype)
saveRDS(Cervical_T, file = "Cervical_T_cells.rds")

###Normalizing the data.
Cervical_T <- NormalizeData(Cervical_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
Cervical_T <- FindVariableFeatures(Cervical_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(Cervical_T)
Cervical_T <- ScaleData(Cervical_T, features = all.genes)
Cervical_T <- ScaleData(Cervical_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
Cervical_T <- RunPCA(Cervical_T, features = VariableFeatures(object = Cervical_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
Cervical_T <- RunHarmony(Cervical_T, group.by.vars = "orig.ident")
Cervical_T <- RunUMAP(Cervical_T, reduction = "harmony", dims = 1:20)
Cervical_T <- RunTSNE(Cervical_T, reduction = "harmony", dims = 1:20)
names(Cervical_T@reductions)

###Cluster the cells.
Cervical_T <- FindNeighbors(Cervical_T, reduction = "harmony", dims = 1:20)
saveRDS(Cervical_T, file = "Cervical_T_cells_before.rds")

#Determine resolution.
test <- Cervical_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_Cervical_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
Cervical_T <- FindClusters(Cervical_T, resolution = 0.1)
saveRDS(Cervical_T, file = "Cervical_T_cells_res0.1.rds")

###classic celltype marker in clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="Cervical_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(Cervical_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="Cervical_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(Cervical_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="Cervical_T_NKT.pdf",width=8, height=5)
DotPlot(Cervical_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="Cervical_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(Cervical_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="Cervical_T_naive.pdf",width=8, height=4)
DotPlot(Cervical_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="Cervical_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(Cervical_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="Cervical_T_Exhaustion.pdf",width=8, height=4)
DotPlot(Cervical_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="Cervical_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(Cervical_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="Cervical_T_NK_marker.pdf",width=8, height=5)
DotPlot(Cervical_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="Cervical_T_cells_T_features_res0.1_annotation.pdf",width=8, height=10)
DotPlot(Cervical_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
CD8_T_cells = c(0)
CD4_T_cells = c(1)
Tregs = c(2)
unknown_T_cells = c(3, 4, 5)
current.cluster.ids <- c(CD8_T_cells,
                         CD4_T_cells,
                         Tregs,
                         unknown_T_cells)
new.cluster.ids <- c(rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("CD4_T_cells",length(CD4_T_cells)),
                     rep("Tregs",length(Tregs)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

Cervical_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Cervical_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Cervical_T$celltype)
table(Cervical_T$seurat_clusters)

###Rename.
Idents(Cervical_T) <- Cervical_T$celltype
Cervical_T <- RenameIdents(Cervical_T, "Tregs" = "Tregs", "CD8_T_cells" = "CD8+ T cells", "CD4_T_cells" = "CD4+ T cells", "unknown_T_cells" = "unknown T cells")
table(Idents(Cervical_T))
Cervical_T$celltype <- Idents(Cervical_T)
table(Cervical_T$celltype)
saveRDS(Cervical_T, file = "Cervical_T_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
Cervical_T <- subset(Cervical_T, subset = celltype %in% c("CD8+ T cells", "CD4+ T cells", "Tregs"))
Cervical_T$seurat_clusters <- as.factor(as.character(Cervical_T$seurat_clusters))
Cervical_T$celltype <- as.factor(as.character(Cervical_T$celltype))
Cervical_T$orig.ident <- as.factor(as.character(Cervical_T$orig.ident))

###Ranking.
Idents(Cervical_T) <- Cervical_T$celltype
Cervical_T <- RenameIdents(Cervical_T, "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(Cervical_T))
Cervical_T$celltype <- Idents(Cervical_T)
table(Cervical_T$celltype)
saveRDS(Cervical_T, file = "Cervical_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6H: left.
pdf(file="Cervical_T_cells_umap_celltype.pdf",width=5.6, height=4)
DimPlot(Cervical_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6H: right.
Idents(Cervical_T) <- 'celltype'
My_levels <- c("Tregs", "CD4+ T cells", "CD8+ T cells")
Idents(Cervical_T) <- factor(Idents(Cervical_T), levels= My_levels)
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3")
pdf(file="Cervical_T_cells_T_features_res0.1_annotation.pdf",width=8, height=5)
DotPlot(Cervical_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
