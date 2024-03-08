###codes for Figure 1G-J.

#1 Figure 1G, H.
#OV.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
###The file GSE165897_UMIcounts_HGSOC.tsv.gz should be extracted to GSE165897_UMIcounts_HGSOC.tsv.
t <- read.delim(file = "GSE165897_UMIcounts_HGSOC.tsv")
t[1:5, 1:5]
rownames(t) <- t[, 1]
t <- t[, -1]
###Move all the parts after the dot to the beginning of the string and connect them with "_" to the original beginning part of the string.
colnames(t) <- gsub("^(.*)\\.(.*)", "\\2_\\1", colnames(t))

###Replace the first "_" in the string with ".", leaving other "_" unchanged.
colnames(t) <- gsub("^([^_]*)_", "\\1.", colnames(t))
saveRDS(t, file = "GSE165897_matrix.rds")

library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)

###Create one Seurat object
OV <- CreateSeuratObject(counts = t, project = "OV")
OV
head(OV@meta.data)

###Annotate the Treatment.
OV$Treatment <- ifelse(test = OV$orig.ident %in% c("EOC372.iPer", "EOC443.iOme1", "EOC540.iOme2", "EOC3.iOme2", "EOC87.iOme1", "EOC136.iOme", "EOC1005.iTum2", "EOC733.iOme", "EOC153.iOme1", "EOC349.iOme1", "EOC227.iOme1"), yes = "NACT", no = "naive")
table(OV$Treatment)

###Annotate Sample.
Idents(OV) <- OV$orig.ident
OV <- RenameIdents(OV, "EOC1005.iTum2" = "EOC1005", "EOC1005.pPer" = "EOC1005", "EOC136.iOme" = "EOC136", "EOC136.pMes1" = "EOC136",
                   "EOC153.iOme1" = "EOC153", "EOC153.pOme" = "EOC153", "EOC227.iOme1" = "EOC227", "EOC227.pOme1" = "EOC227",
                   "EOC3.iOme2" = "EOC3", "EOC3.pPer1" = "EOC3", "EOC349.iOme1" = "EOC349", "EOC349.pPer2" = "EOC349",
                   "EOC372.iPer" = "EOC372", "EOC372.pPer" = "EOC372", "EOC443.iOme1" = "EOC443", "EOC443.pOme" = "EOC443",
                   "EOC540.iOme2" = "EOC540", "EOC540.p2Ome" = "EOC540", "EOC733.iOme" = "EOC733", "EOC733.pPer" = "EOC733",
                   "EOC87.iOme1" = "EOC87", "EOC87.pPer1" = "EOC87")
table(Idents(OV))
OV$Sample <- Idents(OV)
table(OV$Sample)
head(OV@meta.data)
saveRDS(OV, file = "OV_Seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(OV, expression = nFeature_RNA > 200)
selected_f <- rownames(OV)[Matrix::rowSums(OV) > 3]
OV <- subset(OV, features = selected_f, cells = selected_c)
dim(OV)

###QC metrics and filtering.
OV[["percent.mt"]] <- PercentageFeatureSet(OV, pattern = "^MT-")
Idents(OV) <- OV$orig.ident
pdf(file="OV.featureViolin.pdf",width=16,height=6)
VlnPlot(object = OV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="OV.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(OV, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
OV <- subset(OV, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(OV)

###Visualization.
#VlnPlot for Supplementary Figure 1C.
Idents(OV) <- OV$orig.ident
pdf(file="OV.featureViolin_QC.pdf",width=20,height=6)           
VlnPlot(object = OV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
OV <- NormalizeData(OV, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
OV <- FindVariableFeatures(OV, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(OV)
OV <- ScaleData(OV, features = all.genes)
OV <- ScaleData(OV, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
OV <- RunPCA(OV, features = VariableFeatures(object = OV))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
OV <- RunHarmony(OV, group.by.vars = "orig.ident")
OV <- RunUMAP(OV, reduction = "harmony", dims = 1:20)
OV <- RunTSNE(OV, reduction = "harmony", dims = 1:20)
names(OV@reductions)

###Cluster the cells.
OV <- FindNeighbors(OV, reduction = "harmony", dims = 1:20)
saveRDS(OV, file = "OV_before_resolution.rds")

#Determine resolution.
test <- OV
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3B.
library(clustree)
library(patchwork)
pdf(file = "clustree_OV.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
OV <- FindClusters(OV, resolution = 0.2)
saveRDS(OV, file = "OV_resolution_0.2.rds")

###classic celltype marker in clusters.
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
pdf(file="OV_marker.pdf",width=44, height=33)
VlnPlot(OV, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="OV_NK.pdf",width=44, height=14)
VlnPlot(OV, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="OV_Epithelial.pdf",width=44, height=14)
VlnPlot(OV, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="OV_Plasma_cells.pdf",width=44, height=14)
VlnPlot(OV, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="OV_B_marker.pdf",width=44, height=19)
VlnPlot(OV, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="OV_Fibroblasts.pdf",width=44, height=14)
VlnPlot(OV, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="OV_Endothelial.pdf",width=44, height=5)
VlnPlot(OV, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="OV_T_marker.pdf",width=44, height=14)
VlnPlot(OV, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="OV_Mast.pdf",width=44, height=5)
VlnPlot(OV, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="OV_myeloid.pdf",width=44, height=9.5)
VlnPlot(OV, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="OV_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(OV, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="OV_pDC.pdf",width=44, height=20)
VlnPlot(OV, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="OV_neutrophil.pdf",width=44, height=5)
VlnPlot(OV, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="OV_melanocyte.pdf",width=44, height=30)
VlnPlot(OV, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="OV_clear_cell.pdf",width=44, height=10)
VlnPlot(OV, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters
Fibroblasts = c(2)
Myeloid_cells = c(1, 9)
Mast_cells = c(6)
B_cells = c(4)
Plasma_cells = c(7)
T_cells = c(0, 5, 8)
Epithelial_cells = c(3)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Mast_cells,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Epithelial_cells)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells))
)

OV@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(OV@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(OV$Celltype)

###Rename.
Idents(OV) <- OV$Celltype
OV <- RenameIdents(OV, "T_cells" = "T cells", "Plasma_cells" = "Plasma cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells")
table(Idents(OV))
OV$Celltype <- Idents(OV)
table(OV$Celltype)
saveRDS(OV, file = "OV_resolution0.2_annotation.rds")

###Ranking.
Idents(OV) <- OV$Celltype
OV <- RenameIdents(OV, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Epithelial cells" = "Epithelial cells")
table(Idents(OV))
OV$Celltype <- Idents(OV)
table(OV$Celltype)
saveRDS(OV, file = "OV_resolution0.2_annotation2.rds")

###Visualization
#UMAP plot for Figure 1G.
pdf(file="OV_umap_Celltype.pdf",width=5.6, height=4)
DimPlot(OV, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Figure 1H.
Idents(OV) <- OV$Celltype
features = c("KRT19", "EPCAM", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="OV_features.pdf",width=8, height=7.5)
DotPlot(OV, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2B.
pdf(file="OV_umap_orig.ident.pdf",width=7.3, height=4)
DimPlot(OV, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Figure 1I, J.
###Cluster the Fibroblasts.
OV_Fib <- subset(OV, subset = Celltype %in% c("Fibroblasts"))
OV_Fib$seurat_clusters <- as.factor(as.character(OV_Fib$seurat_clusters))
OV_Fib$Celltype <- as.factor(as.character(OV_Fib$Celltype))
OV_Fib$orig.ident <- as.factor(as.character(OV_Fib$orig.ident))
OV_Fib$Treatment <- as.factor(as.character(OV_Fib$Treatment))
saveRDS(OV_Fib, file = "OV_Fib.rds")
Fib <- OV_Fib
table(Fib$Celltype)

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
saveRDS(Fib, file = "OV_Fib_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4B.
library(clustree)
library(patchwork)
pdf(file = "clustree_OV_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
Fib <- FindClusters(Fib, resolution = 0.2)
saveRDS(Fib, file = "OV_Fib_res0.2.rds")

###classic celltype marker in clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="OV_Fib_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="OV_Fib_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="OV_Fib_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="OV_Fib_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5")
pdf(file="OV_Fib_myCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(3)
iCAFs = c(2)
myCAFs = c(4, 6)
eCAFs = c(0)
unknown_CAFs = c(1, 5)
current.cluster.ids <- c(apCAFs,
                         iCAFs,
                         myCAFs,
                         eCAFs,
                         unknown_CAFs)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
                     rep("iCAFs",length(iCAFs)),
                     rep("myCAFs",length(myCAFs)),
                     rep("eCAFs",length(eCAFs)),
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
saveRDS(Fib, file = "OV_Fib_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
Fib <- subset(Fib, subset = celltype %in% c("apCAFs", "iCAFs", "myCAFs", "eCAFs"))
Fib$seurat_clusters <- as.factor(as.character(Fib$seurat_clusters))
Fib$celltype <- as.factor(as.character(Fib$celltype))
Fib$orig.ident <- as.factor(as.character(Fib$orig.ident))
Fib$Treatment <- as.factor(as.character(Fib$Treatment))

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "OV_Fib_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Figure 1I.
pdf(file="OV_Fib_umap_celltype.pdf",width=5.4, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Figure 1J.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "eCAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "CXCL12", "CXCL1", "CCL2", "POSTN", "MMP14", "RGS5", "ACTA2", "HLA-DPB1", "HLA-DRA", "CD74")
pdf(file="OV_Fib_F_features_annotation.pdf",width=8, height=6)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
