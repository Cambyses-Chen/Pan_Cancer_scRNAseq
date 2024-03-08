###codes for Supplementary Figure 6E-F.

#NPC.
#1 Supplementary Figure 6E.
###Cluster the Myeloid cells.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
NPC <- readRDS(file = "NPC_resolution0.2_annotation_sub.rds") #The detailed code for NPC_resolution0.2_annotation_sub.rds is at line 616 of the Supplementary Figure 5A-B.R file.
NPC_M <- subset(NPC, subset = Celltype %in% c("Myeloid cells"))
NPC_M$seurat_clusters <- as.factor(as.character(NPC_M$seurat_clusters))
NPC_M$Celltype <- as.factor(as.character(NPC_M$Celltype))
NPC_M$orig.ident <- as.factor(as.character(NPC_M$orig.ident))
NPC_M$Sample <- as.factor(as.character(NPC_M$Sample))
saveRDS(NPC_M, file = "NPC_Myeloid.rds")
table(NPC_M$Celltype)

###Normalizing the data.
NPC_M <- NormalizeData(NPC_M, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC_M <- FindVariableFeatures(NPC_M, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC_M)
NPC_M <- ScaleData(NPC_M, features = all.genes)
NPC_M <- ScaleData(NPC_M, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC_M <- RunPCA(NPC_M, features = VariableFeatures(object = NPC_M))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC_M <- RunHarmony(NPC_M, group.by.vars = "orig.ident")
NPC_M <- RunUMAP(NPC_M, reduction = "harmony", dims = 1:20)
NPC_M <- RunTSNE(NPC_M, reduction = "harmony", dims = 1:20)
names(NPC_M@reductions)

###Cluster the cells.
NPC_M <- FindNeighbors(NPC_M, reduction = "harmony", dims = 1:20)
saveRDS(NPC_M, file = "NPC_Myeloid_before.rds")

#Determine resolution.
test <- NPC_M
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_NPC_M.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.3.
NPC_M <- FindClusters(NPC_M, resolution = 0.3)
saveRDS(NPC_M, file = "NPC_Myeloid_res0.3.rds")

###classic celltype marker expression across different clusters.
Neutrophil = c("CEACAM8", "FCGR3B", "CSF3R")
pdf(file="NPC_Myeloid_neutrophil.pdf",width=44, height=10)
VlnPlot(NPC_M, features = Neutrophil, pt.size = 0, ncol = 2)
dev.off()

Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1", "MAFB", "C1QA", "ITGAM")
pdf(file="NPC_Myeloid_Monocyte_macrophage.pdf",width=44, height=20)
VlnPlot(NPC_M, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

Macrophage_M2=c("CD14", "CD68", "CD163", "CCL18")
pdf(file="NPC_Myeloid_Macrophage_M2.pdf",width=44, height=9.5)
VlnPlot(NPC_M, features = Macrophage_M2, pt.size = 0, ncol = 2)
dev.off()

pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="NPC_Myeloid_pDC.pdf",width=44, height=20)
VlnPlot(NPC_M, features = pDC, pt.size = 0, ncol = 2)
dev.off()

M1=c("CD80", "CD86", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "MRC1")
pdf(file="NPC_Myeloid_M1.pdf",width=44, height=20)
VlnPlot(NPC_M, features = M1, pt.size = 0, ncol = 2)
dev.off()

Monocyte=c("LYZ", "CD14", "FCGR3A", "FCGR3B", "FCN1", "S100A9")
pdf(file="NPC_Myeloid_Monocyte.pdf",width=44, height=15)
VlnPlot(NPC_M, features = Monocyte, pt.size = 0, ncol = 2)
dev.off()

DC=c("LYZ", "CCR7", "SAMSN1", "FCER1A", "HLA-DQA1")
pdf(file="NPC_Myeloid_DC.pdf",width=44, height=15)
VlnPlot(NPC_M, features = DC, pt.size = 0, ncol = 2)
dev.off()

DC_marker = c("IL12A", "IL6", "TLR7", "TLR9")
pdf(file="NPC_Myeloid_DC_marker2.pdf",width=44, height=9.5)
VlnPlot(NPC_M, features = DC_marker, pt.size = 0, ncol = 2)
dev.off()

cDC1 = c("BTLA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "IDO1", "HLA-DQA1")
pdf(file="NPC_Myeloid_cDC1_marker.pdf",width=44, height=30)
VlnPlot(NPC_M, features = cDC1, pt.size = 0, ncol = 2)
dev.off()

cDC2 = c("CD1C", "CD1E", "CD1A", "ITGAX", "ITGAM", "FCER1A", "SIRPA", "CLEC10A", "CD163", "CD2", "CX3CR1", "CD14", "NOTCH2", "BATF3", "ID2", "IRF8", "ZBTB46", "HLA-DQA1")
pdf(file="NPC_Myeloid_cDC2_marker.pdf",width=44, height=45)
VlnPlot(NPC_M, features = cDC2, pt.size = 0, ncol = 2)
dev.off()

moDC = c("CD14", "CD1A", "CD1C", "CD209", "FCER1A", "ITGAM", "MRC1", "SIRPA", "FCGR3A", "S100A8", "S100A9")
pdf(file="NPC_Myeloid_moDC_marker.pdf",width=44, height=30)
VlnPlot(NPC_M, features = moDC, pt.size = 0, ncol = 2)
dev.off()

LC = c("CD1A", "CD207", "ID2")
pdf(file="NPC_Myeloid_LC_marker.pdf",width=44, height=10)
VlnPlot(NPC_M, features = LC, pt.size = 0, ncol = 2)
dev.off()

preDC = c("CD5", "AXL", "SIGLEC6")
pdf(file="NPC_Myeloid_preDC_marker.pdf",width=44, height=10)
VlnPlot(NPC_M, features = preDC, pt.size = 0, ncol = 2)
dev.off()

mast=c("KIT", "TPSAB1")
pdf(file="NPC_Myeloid_mast.pdf",width=44, height=5)
VlnPlot(NPC_M, features = mast, pt.size = 0, ncol = 2)
dev.off()

M_features = c("FCGR3B", "KATNBL1", "CD163", "CCR7", "LAMP3", "CD14", "CD68", "CLEC10A", "CD1C", "CLEC9A", "CADM1", "SPP1", "C1QC", "C1QA", "APOE", "MRC1", "VCAN", "THBS1", "FCGR3A", "S100A8", "S100A9", "FCN1", "MKI67", "LYVE1", "PLTP", "SEPP1", "INHBA", "CCL4", "IL1RN", "NLRP3", "EREG", "IL1B", "LILRA4", "PTGDS", "PLAC8", "FCER1A")
pdf(file="NPC_Myeloid_features.pdf",width=10, height=14)
DotPlot(NPC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
cDC1 = c(9)
cDC2 = c(2)
LAMP3_DCs = c(6)
Langerhans_cells = c(7)
VCAN_Macrophages = c(0)
C1QC_Macrophages = c(1, 4)
unknown_Myeloid = c(3, 5, 8, 10, 11)
current.cluster.ids <- c(cDC1,
                         cDC2,
                         LAMP3_DCs,
                         Langerhans_cells,
                         VCAN_Macrophages,
                         C1QC_Macrophages,
                         unknown_Myeloid)
new.cluster.ids <- c(rep("cDC1",length(cDC1)),
                     rep("cDC2",length(cDC2)),
                     rep("LAMP3_DCs",length(LAMP3_DCs)),
                     rep("Langerhans_cells",length(Langerhans_cells)),
                     rep("VCAN_Macrophages",length(VCAN_Macrophages)),
                     rep("C1QC_Macrophages",length(C1QC_Macrophages)),
                     rep("unknown_Myeloid",length(unknown_Myeloid))
                     
)

NPC_M@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(NPC_M@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(NPC_M$celltype)
table(NPC_M$seurat_clusters)

###Rename.
Idents(NPC_M) <- NPC_M$celltype
NPC_M <- RenameIdents(NPC_M, "LAMP3_DCs" = "LAMP3+ DCs", "Langerhans_cells" = "Langerhans cells", "C1QC_Macrophages" = "C1QC+ Macrophages", "VCAN_Macrophages" = "VCAN+ Macrophages", "unknown_Myeloid" = "unknown Myeloid cells")
table(Idents(NPC_M))
NPC_M$celltype <- Idents(NPC_M)
table(NPC_M$celltype)
saveRDS(NPC_M, file = "NPC_Myeloid_res0.3_annotation.rds")

###Exclude unknown cell subpopulations.
NPC_M <- subset(NPC_M, subset = celltype %in% c("LAMP3+ DCs", "C1QC+ Macrophages", "VCAN+ Macrophages", "cDC1", "cDC2", "Langerhans cells"))
NPC_M$seurat_clusters <- as.factor(as.character(NPC_M$seurat_clusters))
NPC_M$celltype <- as.factor(as.character(NPC_M$celltype))
NPC_M$orig.ident <- as.factor(as.character(NPC_M$orig.ident))

###Ranking.
Idents(NPC_M) <- NPC_M$celltype
NPC_M <- RenameIdents(NPC_M, "cDC1" = "cDC1", "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "Langerhans cells" = "Langerhans cells", "C1QC+ Macrophages" = "C1QC+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(NPC_M))
NPC_M$celltype <- Idents(NPC_M)
table(NPC_M$celltype)
saveRDS(NPC_M, file = "NPC_Myeloid_res0.3_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6E: left.
pdf(file="NPC_Myeloid_umap_celltype.pdf",width=6.2, height=4)
DimPlot(NPC_M, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6E: right.
Idents(NPC_M) <- 'celltype'
My_levels <- c( "cDC1", "cDC2", "LAMP3+ DCs", "Langerhans cells", "C1QC+ Macrophages", "VCAN+ Macrophages")
Idents(NPC_M) <- factor(Idents(NPC_M), levels= My_levels)
M_features = c("FCGR3A", "CD163", "CD14", "CD68", "VCAN", "C1QC", "CD1A", "CD207", "CCR7", "LAMP3", "CLEC10A", "CD1C", "CLEC9A", "CADM1")
pdf(file="NPC_Myeloid_Myeloid_features.pdf",width=8, height=8)
DotPlot(NPC_M, features = M_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#2 Supplementary Figure 6F.
###Cluster the T cells.
NPC_T <- subset(NPC, subset = Celltype %in% c("T cells"))
NPC_T$seurat_clusters <- as.factor(as.character(NPC_T$seurat_clusters))
NPC_T$Celltype <- as.factor(as.character(NPC_T$Celltype))
NPC_T$orig.ident <- as.factor(as.character(NPC_T$orig.ident))
NPC_T$Sample <- as.factor(as.character(NPC_T$Sample))
saveRDS(NPC_T, file = "NPC_T.rds")

###Normalizing the data.
NPC_T <- NormalizeData(NPC_T, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC_T <- FindVariableFeatures(NPC_T, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC_T)
NPC_T <- ScaleData(NPC_T, features = all.genes)
NPC_T <- ScaleData(NPC_T, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC_T <- RunPCA(NPC_T, features = VariableFeatures(object = NPC_T))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC_T <- RunHarmony(NPC_T, group.by.vars = "orig.ident")
NPC_T <- RunUMAP(NPC_T, reduction = "harmony", dims = 1:20)
NPC_T <- RunTSNE(NPC_T, reduction = "harmony", dims = 1:20)
names(NPC_T@reductions)

###Cluster the cells.
NPC_T <- FindNeighbors(NPC_T, reduction = "harmony", dims = 1:20)
saveRDS(NPC_T, file = "NPC_T_cells_before.rds")

#Determine resolution.
test <- NPC_T
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}
library(clustree)
library(patchwork)
pdf(file = "clustree_NPC_T.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
NPC_T <- FindClusters(NPC_T, resolution = 0.2)
saveRDS(NPC_T, file = "NPC_T_cells_res0.2.rds")

###classic celltype marker expression across different clusters.
CD8_T_cells=c("CD3D", "CD8A", "CD8B")
pdf(file="NPC_T_CD8_T_cells.pdf",width=8, height=4)
DotPlot(NPC_T, features = CD8_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_T_cells=c("CD3D", "CD4", "LDHB", "IL7R")
pdf(file="NPC_T_CD4_T_cells.pdf",width=8, height=4)
DotPlot(NPC_T, features = CD4_T_cells, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NKT=c("CD3D", "NKG7", "KLRB1", "KLRD1", "FCGR3A", "FCGR3B", "NCAM1")
pdf(file="NPC_T_NKT.pdf",width=8, height=5)
DotPlot(NPC_T, features = NKT, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

CD4_Treg=c("IKZF2", "IL2RA", "FOXP3", "CD4")
pdf(file="NPC_T_CD4_Treg.pdf",width=8, height=4)
DotPlot(NPC_T, features = CD4_Treg, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

naive = c("CCR7", "LEF1", "SELL", "TCF7", "IL7R", "CD4")
pdf(file="NPC_T_naive.pdf",width=8, height=4)
DotPlot(NPC_T, features = naive, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Cytotoxic=c("IFNG", "PRF1", "GZMK", "GNLY", "GZMA", "NKG7", "CD8A", "CD8B")
pdf(file="NPC_T_Cytotoxic.pdf",width=8, height=5)
DotPlot(NPC_T, features = Cytotoxic, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Exhaustion=c("CTLA4", "HAVCR2", "LAG3", "TIGIT", "PDCD1")
pdf(file="NPC_T_Exhaustion.pdf",width=8, height=4)
DotPlot(NPC_T, features = Exhaustion, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

Co_simulatory=c("TNFRSF9", "ICOS", "TNFRSF14", "CD28")
pdf(file="NPC_T_Co_simulatory.pdf",width=8, height=4)
DotPlot(NPC_T, features = Co_simulatory, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="NPC_T_NK_marker.pdf",width=8, height=5)
DotPlot(NPC_T, features = NK_marker, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

T_features = c("CD8A", "CD8B", "CD3D", "CD3E", "CD4", "IFNG", "GZMK", "GNLY", "GZMA", "NKG7", "PRF1", "CCR7", "SELL", "IL7R", "TCF7", "IL17A", "IL22", "CXCR5", "BCL6", "CXCR3", "CCR5", "IL4", "IL2RA", "FOXP3")
pdf(file="NPC_T_cells_T_features.pdf",width=8, height=10)
DotPlot(NPC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
Tregs = c(3)
CD8_T_cells = c(1)
NK_cells = c(6)
CD4_T_cells = c(0)
unknown_T_cells = c(2, 4, 5, 7)
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

NPC_T@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(NPC_T@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(NPC_T$seurat_clusters)
table(NPC_T$celltype)

###Rename.
Idents(NPC_T) <- NPC_T$celltype
NPC_T <- RenameIdents(NPC_T, "CD8_T_cells" = "CD8+ T cells", "NK_cells" = "NK cells", "CD4_T_cells" = "CD4+ T cells", "unknown_T_cells" = "unknown T cells")
table(Idents(NPC_T))
NPC_T$celltype <- Idents(NPC_T)
table(NPC_T$celltype)
saveRDS(NPC_T, file = "NPC_T_cells_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
NPC_T <- subset(NPC_T, subset = celltype %in% c("CD8+ T cells", "NK cells", "Tregs", "CD4+ T cells"))
NPC_T$seurat_clusters <- as.factor(as.character(NPC_T$seurat_clusters))
NPC_T$celltype <- as.factor(as.character(NPC_T$celltype))

###Ranking.
Idents(NPC_T) <- NPC_T$celltype
NPC_T <- RenameIdents(NPC_T, "NK cells" = "NK cells", "Tregs" = "Tregs", "CD4+ T cells" = "CD4+ T cells", "CD8+ T cells" = "CD8+ T cells")
table(Idents(NPC_T))
NPC_T$celltype <- Idents(NPC_T)
table(NPC_T$celltype)
saveRDS(NPC_T, file = "NPC_T_cells_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Supplementary Figure 6F: left.
pdf(file="NPC_T_cells_umap_celltype.pdf",width=5.6, height=4)
DimPlot(NPC_T, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 6F: right.
Idents(NPC_T) <- NPC_T$celltype
T_features = c("CD8B", "CD8A", "CD4", "IL2RA", "FOXP3", "FGFBP2", "KLRD1")
pdf(file="NPC_T_cells_T_features_res0.2_annotation.pdf",width=8, height=4.5)
DotPlot(NPC_T, features = T_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
