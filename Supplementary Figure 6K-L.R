###codes for Supplementary Figure 6K-L.

#BRCA.
#1 Supplementary Figure 6K.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
BC <- readRDS(file = "BC_GSE176078_resolution0.1_annotation_sub.rds") #The detailed code for BC_GSE176078_resolution0.1_annotation_sub.rds is at line 308 of the Supplementary Figure 5G-H.R file.
BC_M <- subset(BC, subset = Celltype %in% c("Myeloid cells"))
table(BC_M$Celltype)
table(BC_M$orig.ident)
BC_M$seurat_clusters <- as.factor(as.character(BC_M$seurat_clusters))
BC_M$Celltype <- as.factor(as.character(BC_M$Celltype))
saveRDS(BC_M, file = "BC_GSE176078_Myeloid.rds")

###Normalizing the data.
BC_M <- NormalizeData(BC_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
BC_M <- FindVariableFeatures(BC_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(BC_M)
BC_M <- ScaleData(BC_M, features = all.genes)
BC_M <- ScaleData(BC_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
BC_M <- RunPCA(BC_M, features = VariableFeatures(object = BC_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
BC_M <- RunHarmony(BC_M, group.by.vars = "orig.ident")
BC_M <- RunUMAP(BC_M, reduction = "harmony", dims = 1:20)
BC_M <- RunTSNE(BC_M, reduction = "harmony", dims = 1:20)
names(BC_M@reductions)

###Cluster the cells.
BC_M <- FindNeighbors(BC_M, reduction = "harmony", dims = 1:20)
saveRDS(BC_M, file = "BC_GSE176078_Myeloid_before.rds")

#Determine resolution.
test <- BC_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_BC_Myeloid.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
BC_M <- FindClusters(BC_M, resolution = 0.2)
saveRDS(BC_M, file = "BC_GSE176078_Myeloid_res0.2.rds")

###classic celltype marker in clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="BC_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(BC_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="BC_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(BC_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="BC_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(BC_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="BC_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(BC_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="BC_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(BC_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="BC_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(BC_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="BC_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(BC_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="BC_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(BC_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="BC_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(BC_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="BC_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(BC_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="BC_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(BC_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="BC_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(BC_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="BC_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(BC_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="BC_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(BC_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="BC_Myeloid_features.pdf",width=10, height=14)
DotPlot(BC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(6)
cDC2 = c(2)
LAMP3_DCs = c(7)
C1QC_Macrophages = c(0, 1)
SPP1_Macrophages = c(5)
VCAN_Macrophages = c(4)
unknown_Myeloid = c(3, 8)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         LAMP3_DCs,
                         C1QC_Macrophages,
                         SPP1_Macrophages,
                         VCAN_Macrophages,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("SPP1_Macrophages",length(SPP1_Macrophages)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
                     
)

BC_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(BC_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(BC_M$celltype)
table(BC_M$seurat_clusters)

###Rename.
Idents(BC_M) <- BC_M$celltype
BC_M <- RenameIdents(BC_M, "LAMP3_DCs" = "LAMP3+ DCs", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages", "SPP1_Macrophages" = "SPP1+ Macrophages", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(BC_M))
BC_M$celltype <- Idents(BC_M)
table(BC_M$celltype)
saveRDS(BC_M, file = "BC_GSE176078_Myeloid_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
BC_M <- subset(BC_M, subset = celltype %in% c("LAMP3+ DCs", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages", "cDC1", "cDC2"))
BC_M$seurat_clusters <- as.factor(as.character(BC_M$seurat_clusters))
BC_M$celltype <- as.factor(as.character(BC_M$celltype))

###Ranking.
Idents(BC_M) <- BC_M$celltype
BC_M <- RenameIdents(BC_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(BC_M))
BC_M$celltype <- Idents(BC_M)
table(BC_M$celltype)
saveRDS(BC_M, file = "BC_GSE176078_Myeloid_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6K: left.
pdf(file="BC_GSE176078_Myeloid_umap_celltype.pdf",width=6, height=4)
DimPlot(BC_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6K: right.
Idents(BC_M) <- 'celltype'
My_levels <- c( "cDC1", "cDC2", "LAMP3+ DCs", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages")
Idents(BC_M) <- factor(Idents(BC_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "SPP1", "C1QC", "CCR7", "LAMP3", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="BC_GSE176078_Myeloid_M_features_res0.2.pdf",width=8, height=7.5)
DotPlot(BC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6L.
###Cluster the T cells.
BC_T <- subset(BC, subset = Celltype %in% c("T cells"))
table(BC_T$Celltype)
BC_T$seurat_clusters <- as.factor(as.character(BC_T$seurat_clusters))
BC_T$Celltype <- as.factor(as.character(BC_T$Celltype))
saveRDS(BC_T, file = "BC_GSE176078_T_cells.rds")

###Normalizing the data.
BC_T <- NormalizeData(BC_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
BC_T <- FindVariableFeatures(BC_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(BC_T)
BC_T <- ScaleData(BC_T, features = all.genes)
BC_T <- ScaleData(BC_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
BC_T <- RunPCA(BC_T, features = VariableFeatures(object = BC_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
BC_T <- RunHarmony(BC_T, group.by.vars = "orig.ident")
BC_T <- RunUMAP(BC_T, reduction = "harmony", dims = 1:20)
BC_T <- RunTSNE(BC_T, reduction = "harmony", dims = 1:20)
names(BC_T@reductions)

###Cluster the cells.
BC_T <- FindNeighbors(BC_T, reduction = "harmony", dims = 1:20)
saveRDS(BC_T, file = "BC_GSE176078_T_cells_before.rds")

#Determine resolution.
test <- BC_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_BC_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
BC_T <- FindClusters(BC_T, resolution = 0.1)
saveRDS(BC_T, file = "BC_GSE176078_T_cells_res0.1.rds")

###classic celltype marker in clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="BC_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(BC_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="BC_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(BC_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="BC_T_NKT.pdf",width=8, height=5)
DotPlot(BC_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="BC_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(BC_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="BC_T_naive.pdf",width=8, height=4)
DotPlot(BC_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="BC_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(BC_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="BC_T_Exhaustion.pdf",width=8, height=4)
DotPlot(BC_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="BC_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(BC_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="BC_T_NK_marker.pdf",width=8, height=5)
DotPlot(BC_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="BC_T_cells_T_features_res0.1_annotation.pdf",width=8, height=10)
DotPlot(BC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(2)
CD8_T_cells = c(0)
NK_cells = c(4)
CD4_T_cells = c(1)
unknown_T_cells = c(3, 5)
current.cluster.ids <- c(Tregs,
                         CD8_T_cells,
                         NK_cells,
                         CD4_T_cells,
                         unknown_T_cells)
new.cluster.ids <- c(rep("Tregs",length(Tregs)),
                     rep("CD8_T_cells",length(CD8_T_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("CD4_T_cells",length(CD4_T_cells)),
                     rep("unknown_T_cells",length(unknown_T_cells))
)

BC_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(BC_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(BC_T$seurat_clusters)
table(BC_T$celltype)

###Rename.
Idents(BC_T) <- BC_T$celltype
BC_T <- RenameIdents(BC_T, "CD8_T_cells" = "CD8+ T cells", "CD4_T_cells" = "CD4+ T cells", "NK_cells" = "NK cells", "unknown_T_cells" = "unknown T cells")
table(Idents(BC_T))
BC_T$celltype <- Idents(BC_T)
table(BC_T$celltype)
saveRDS(BC_T, file = "BC_GSE176078_T_cells_res0.1_annotation.rds")

###Exclude unknown cell subpopulations.
BC_T <- subset(BC_T, subset = celltype %in% c("CD8+ T cells", "CD4+ T cells", "NK cells", "Tregs"))
BC_T$seurat_clusters <- as.factor(as.character(BC_T$seurat_clusters))
BC_T$celltype <- as.factor(as.character(BC_T$celltype))
BC_T$orig.ident <- as.factor(as.character(BC_T$orig.ident))

###Ranking.
Idents(BC_T) <- BC_T$celltype
BC_T <- RenameIdents(BC_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(BC_T))
BC_T$celltype <- Idents(BC_T)
table(BC_T$celltype)
saveRDS(BC_T, file = "BC_GSE176078_T_cells_res0.1_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6L: left.
pdf(file="BC_GSE176078_T_cells_umap_celltype.pdf",width=5.3, height=4)
DimPlot(BC_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6L: right.
Idents(BC_T) <- BC_T$celltype
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="BC_GSE176078_T_cells_T_features_res0.1_annotation.pdf",width=8, height=4.5)
DotPlot(BC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
