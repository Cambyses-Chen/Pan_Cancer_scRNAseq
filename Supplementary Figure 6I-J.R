###codes for Supplementary Figure 6I-J.

#CRC.
#1 Supplementary Figure 6I.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
CRC <- readRDS(file = "CRC_HRA000979_resolution0.1_annotation_sub.rds") #The detailed code for CRC_HRA000979_resolution0.1_annotation_sub.rds is at line 269 of the Supplementary Figure 5E-F.R file.
CRC_M <- subset(CRC, subset = Celltype %in% c("Myeloid cells"))
CRC_M$seurat_clusters <- as.factor(as.character(CRC_M$seurat_clusters))
CRC_M$Celltype <- as.factor(as.character(CRC_M$Celltype))
CRC_M$orig.ident <- as.factor(as.character(CRC_M$orig.ident))
CRC_M$Sample <- as.factor(as.character(CRC_M$Sample))
table(CRC_M$Celltype)
saveRDS(CRC_M, file = "CRC_HRA000979_Myeloid.rds")

###Normalizing the data.
CRC_M <- NormalizeData(CRC_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CRC_M <- FindVariableFeatures(CRC_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CRC_M)
CRC_M <- ScaleData(CRC_M, features = all.genes)
CRC_M <- ScaleData(CRC_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CRC_M <- RunPCA(CRC_M, features = VariableFeatures(object = CRC_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CRC_M <- RunHarmony(CRC_M, group.by.vars = "orig.ident")
CRC_M <- RunUMAP(CRC_M, reduction = "harmony", dims = 1:20)
CRC_M <- RunTSNE(CRC_M, reduction = "harmony", dims = 1:20)
names(CRC_M@reductions)

###Cluster the cells.
CRC_M <- FindNeighbors(CRC_M, reduction = "harmony", dims = 1:20)
saveRDS(CRC_M, file = "CRC_HRA000979_Myeloid_before.rds")

#Determine resolution.
test <- CRC_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_CRC_Mye.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
CRC_M <- FindClusters(CRC_M, resolution = 0.2)
saveRDS(CRC_M, file = "CRC_HRA000979_Myeyeloid_res0.2.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="CRC_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(CRC_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="CRC_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(CRC_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="CRC_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(CRC_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="CRC_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(CRC_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="CRC_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(CRC_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="CRC_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(CRC_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="CRC_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(CRC_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="CRC_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(CRC_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="CRC_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(CRC_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="CRC_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(CRC_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="CRC_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(CRC_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="CRC_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(CRC_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="CRC_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(CRC_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="CRC_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(CRC_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="CRC_Myeloid_features.pdf",width=10, height=14)
DotPlot(CRC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC2 = c(4)
VCAN_Macrophages = c(1)
C1QC_Macrophages = c(0, 3)
Neutrophils = c(2)
LAMP3_DCs = c(5)
unknown_Myeloid = c(6)
current.cluster.ids <- c(cDC2,
                         VCAN_Macrophages,
                         C1QC_Macrophages,
                         Neutrophils,
                         LAMP3_DCs,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC2",length(cDC2)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("Neutrophils",length(Neutrophils)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
)

CRC_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(CRC_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CRC_M$celltype)
table(CRC_M$seurat_clusters)

###Rename.
Idents(CRC_M) <- CRC_M$celltype
CRC_M <- RenameIdents(CRC_M, "LAMP3_DCs" = "LAMP3+ DCs", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(CRC_M))
CRC_M$celltype <- Idents(CRC_M)
table(CRC_M$celltype)
saveRDS(CRC_M, file = "CRC_HRA000979_Myeloid_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
CRC_M <- subset(CRC_M, subset = celltype %in% c("LAMP3+ DCs", "C1QC+ Macrophages", "VCAN+ Macrophages", "cDC2", "Neutrophils"))
CRC_M$seurat_clusters <- as.factor(as.character(CRC_M$seurat_clusters))
CRC_M$celltype <- as.factor(as.character(CRC_M$celltype))
CRC_M$orig.ident <- as.factor(as.character(CRC_M$orig.ident))

###Ranking.
Idents(CRC_M) <- CRC_M$celltype
CRC_M <- RenameIdents(CRC_M, "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "Neutrophils" = "Neutrophils", "C1QC+ Macrophages" = "C1QC+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(CRC_M))
CRC_M$celltype <- Idents(CRC_M)
table(CRC_M$celltype)
saveRDS(CRC_M, file = "CRC_HRA000979_Myeloid_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6I: left.
pdf(file="CRC_HRA000979_Myeloid_umap_celltype.pdf",width=6, height=4)
DimPlot(CRC_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv() 
dev.off()

#Bubble heatmap for Supplementary Figure 6I: right.
Idents(CRC_M) <- 'celltype'
My_levels <- c( "cDC2", "LAMP3+ DCs", "Neutrophils", "C1QC+ Macrophages", "VCAN+ Macrophages")
Idents(CRC_M) <- factor(Idents(CRC_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "C1QC", "KATNBL1", "FCGR3B", "CCR7", "LAMP3", "CLEC10A", "CD1C")
pdf(file="CRC_HRA000979_Myeloid_M_features_res0.2.pdf",width=8, height=7.5)
DotPlot(CRC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6J.
###Cluster the T cells.
CRC_T <- subset(CRC, subset = Celltype %in% c("T cells"))
CRC_T$seurat_clusters <- as.factor(as.character(CRC_T$seurat_clusters))
CRC_T$Celltype <- as.factor(as.character(CRC_T$Celltype))
CRC_T$orig.ident <- as.factor(as.character(CRC_T$orig.ident))
CRC_T$Sample <- as.factor(as.character(CRC_T$Sample))
table(CRC_T$Celltype)
saveRDS(CRC_T, file = "CRC_HRA000979_T_cells.rds")

###Normalizing the data.
CRC_T <- NormalizeData(CRC_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CRC_T <- FindVariableFeatures(CRC_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CRC_T)
CRC_T <- ScaleData(CRC_T, features = all.genes)
CRC_T <- ScaleData(CRC_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CRC_T <- RunPCA(CRC_T, features = VariableFeatures(object = CRC_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CRC_T <- RunHarmony(CRC_T, group.by.vars = "orig.ident")
CRC_T <- RunUMAP(CRC_T, reduction = "harmony", dims = 1:20)
CRC_T <- RunTSNE(CRC_T, reduction = "harmony", dims = 1:20)
names(CRC_T@reductions)

###Cluster the cells.
CRC_T <- FindNeighbors(CRC_T, reduction = "harmony", dims = 1:20)
saveRDS(CRC_T, file = "CRC_HRA000979_T_cells_before.rds")

#Determine resolution.
test <- CRC_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_CRC_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
CRC_T <- FindClusters(CRC_T, resolution = 0.1)
saveRDS(CRC_T, file = "CRC_HRA000979_T_cells_res0.1.rds")

###classic celltype marker expression across different clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="CRC_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(CRC_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="CRC_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(CRC_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="CRC_T_NKT.pdf",width=8, height=5)
DotPlot(CRC_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="CRC_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(CRC_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4", "CD44", "CD69", "IL2RA")
pdf(file="CRC_T_naive.pdf",width=8, height=4)
DotPlot(CRC_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="CRC_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(CRC_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="CRC_T_Exhaustion.pdf",width=8, height=4)
DotPlot(CRC_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="CRC_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(CRC_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="CRC_T_NK_marker.pdf",width=8, height=5)
DotPlot(CRC_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "SELL", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="CRC_T_cells_T_features.pdf",width=8, height=10)
DotPlot(CRC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(2)
CD8_T_cells = c(1)
CD4_T_cells = c(0)
NK_cells = c(3)
unknown_T_cells = c(4)
current.cluster.ids <- c(Tregs,
                         CD8_T_cells,
                         CD4_T_cells,
                         NK_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("Tregs",length(Tregs)),
                     rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("CD4_T_cells",length(CD4_T_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

CRC_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(CRC_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CRC_T$seurat_clusters)
table(CRC_T$celltype)

###Rename.
Idents(CRC_T) <- CRC_T$celltype
CRC_T <- RenameIdents(CRC_T, "CD8_T_cells" = "CD8+ T cells", "CD4_T_cells" = "CD4+ T cells", "NK_cells" = "NK cells", "unknown_T_cells" = "unknown T cells")
table(Idents(CRC_T))
CRC_T$celltype <- Idents(CRC_T)
table(CRC_T$celltype)
saveRDS(CRC_T, file = "CRC_HRA000979_T_cells_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
CRC_T <- subset(CRC_T, subset = celltype %in% c("CD8+ T cells", "CD4+ T cells",  "NK cells", "Tregs"))
CRC_T$seurat_clusters <- as.factor(as.character(CRC_T$seurat_clusters))
CRC_T$celltype <- as.factor(as.character(CRC_T$celltype))
CRC_T$orig.ident <- as.factor(as.character(CRC_T$orig.ident))

###Ranking.
Idents(CRC_T) <- CRC_T$celltype
CRC_T <- RenameIdents(CRC_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(CRC_T))
CRC_T$celltype <- Idents(CRC_T)
table(CRC_T$celltype)
saveRDS(CRC_T, file = "CRC_HRA000979_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6J: left.
pdf(file="CRC_HRA000979_T_cells_umap_celltype.pdf",width=5.5, height=4)
DimPlot(CRC_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6J: right.
Idents(CRC_T) <- 'celltype'
My_levels <- c("NK cells", "Tregs", "CD4+ T cells", "CD8+ T cells")
Idents(CRC_T) <- factor(Idents(CRC_T), levels= My_levels)
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="CRC_T_cells_T_features_res0.1_annotation.pdf",width=8, height=4.5)
DotPlot(CRC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
