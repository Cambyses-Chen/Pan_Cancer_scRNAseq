###codes for Figure 1C-F.

#1 Figure 1C, D.
#HNSCC.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
###The file GSE181919_UMI_counts.zip should be extracted to GSE181919_UMI_counts.txt. The file GSE181919_Barcode_metadata.txt.gz should be extracted to GSE181919_Barcode_metadata.txt.
t <- read.delim(file = "GSE181919_UMI_counts.txt")
b <- read.delim(file = "GSE181919_Barcode_metadata.txt")
b$barcode <- gsub("-", ".", rownames(b))
identical(colnames(t), b$barcode)
b$barcode2 <- paste0(b$sample.id, "_", b$barcode)
colnames(t) <- b$barcode2
###Exclude samples of Leukoplakia, starting with LP.
t2 <- t[, !grepl("^LP", names(t))]

###Create one Seurat object
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
HNSCC <- CreateSeuratObject(counts = t2, project = "HNSCC")
table(HNSCC$orig.ident)
HNSCC
head(colnames(HNSCC))
tail(colnames(HNSCC))

library(dplyr)
#Comment sample information.
HNSCC@meta.data <- HNSCC@meta.data %>%
  mutate(
    Sample = case_when(
      orig.ident %in% c("C04", "C06", "C07", "C08", "C09", "C15", "C21", "C22", "C26", "C30", "C31", "C38", "C43", "C46", "C51", "C57", "C59", "C60", "C84", "C86") ~ "Primary",
      orig.ident %in% c("LN22", "LN38", "LN46", "LN59") ~ "Metastatic",
      TRUE ~ "Normal"  #If all conditions are not met, set it to "Normal"
    )
  )
table(HNSCC$Sample)
saveRDS(HNSCC, file = "HNSCC_Seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(HNSCC, expression = nFeature_RNA > 200)
selected_f <- rownames(HNSCC)[Matrix::rowSums(HNSCC) > 3]
HNSCC <- subset(HNSCC, features = selected_f, cells = selected_c)
dim(HNSCC)

###QC metrics and filtering.
HNSCC[["percent.mt"]] <- PercentageFeatureSet(HNSCC, pattern = "^MT-")
Idents(HNSCC) <- HNSCC$orig.ident
pdf(file="HNSCC.featureViolin.pdf",width=28,height=6)           
VlnPlot(object = HNSCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="HNSCC.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(HNSCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HNSCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
HNSCC <- subset(HNSCC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(HNSCC)

###Visualization.
#VlnPlot for Supplementary Figure 1A.
Idents(HNSCC) <- HNSCC$orig.ident
pdf(file="HNSCC.featureViolin_QC.pdf",width=24,height=6)           
VlnPlot(object = HNSCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
HNSCC <- NormalizeData(HNSCC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
HNSCC <- FindVariableFeatures(HNSCC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(HNSCC)
HNSCC <- ScaleData(HNSCC, features = all.genes)
HNSCC <- ScaleData(HNSCC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
HNSCC <- RunPCA(HNSCC, features = VariableFeatures(object = HNSCC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
HNSCC$orig.ident <- as.factor(as.character(HNSCC$orig.ident))
HNSCC <- RunHarmony(HNSCC, group.by.vars = "orig.ident")
HNSCC <- RunUMAP(HNSCC, reduction = "harmony", dims = 1:20)
HNSCC <- RunTSNE(HNSCC, reduction = "harmony", dims = 1:20)
names(HNSCC@reductions)

###Cluster the cells.
HNSCC <- FindNeighbors(HNSCC, reduction = "harmony", dims = 1:20)
saveRDS(HNSCC, file = "HNSCC_before_resolution.rds")

#Determine resolution.
test <- HNSCC
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3A.
library(clustree)
library(patchwork)
pdf(file = "clustree_HNSCC.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.3.
HNSCC <- FindClusters(HNSCC, resolution = 0.3)
saveRDS(HNSCC, file = "HNSCC_resolution_0.3.rds")

###classic celltype marker expression across different clusters.
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
pdf(file="HNSCC_marker.pdf",width=44, height=33)
VlnPlot(HNSCC, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="HNSCC_NK.pdf",width=44, height=14)
VlnPlot(HNSCC, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="HNSCC_Epithelial.pdf",width=44, height=14)
VlnPlot(HNSCC, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="HNSCC_Plasma_cells.pdf",width=44, height=14)
VlnPlot(HNSCC, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="HNSCC_B_marker.pdf",width=44, height=19)
VlnPlot(HNSCC, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="HNSCC_Fibroblasts.pdf",width=44, height=14)
VlnPlot(HNSCC, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="HNSCC_Endothelial.pdf",width=44, height=5)
VlnPlot(HNSCC, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="HNSCC_T_marker.pdf",width=44, height=14)
VlnPlot(HNSCC, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="HNSCC_Mast.pdf",width=44, height=5)
VlnPlot(HNSCC, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="HNSCC_myeloid.pdf",width=44, height=9.5)
VlnPlot(HNSCC, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="HNSCC_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(HNSCC, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="HNSCC_pDC.pdf",width=44, height=20)
VlnPlot(HNSCC, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="HNSCC_neutrophil.pdf",width=44, height=5)
VlnPlot(HNSCC, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="HNSCC_melanocyte.pdf",width=44, height=30)
VlnPlot(HNSCC, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="HNSCC_clear_cell.pdf",width=44, height=10)
VlnPlot(HNSCC, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(0, 10, 11)
Myeloid_cells = c(3, 12)
Mast_cells = c(14)
pDCs = c(9)
B_cells = c(4)
Plasma_cells = c(6)
T_cells = c(1, 2, 8)
Epithelial_cells = c(5, 15, 16)
Endothelial_cells = c(7)
unknown = c(13)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Mast_cells,
                         pDCs,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Epithelial_cells,
                         Endothelial_cells,
                         unknown)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("pDCs",length(pDCs)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("unknown",length(unknown))
)

HNSCC@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(HNSCC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(HNSCC$Celltype)
table(HNSCC$seurat_clusters)

###Rename.
Idents(HNSCC) <- HNSCC$Celltype
HNSCC <- RenameIdents(HNSCC, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells", "Plasma_cells" = "Plasma cells")
table(Idents(HNSCC))
HNSCC$Celltype <- Idents(HNSCC)
table(HNSCC$Celltype)
saveRDS(HNSCC, file = "HNSCC_resolution0.3_annotation.rds")

###Exclude unknown cell subpopulations.
HNSCC <- subset(HNSCC, subset = Celltype %in% c("T cells", "Endothelial cells", "B cells", "Epithelial cells", "Myeloid cells", "Mast cells", "Plasma cells", "Fibroblasts", "pDCs"))
HNSCC$seurat_clusters <- as.factor(as.character(HNSCC$seurat_clusters))
HNSCC$Celltype <- as.factor(as.character(HNSCC$Celltype))
HNSCC$orig.ident <- as.factor(as.character(HNSCC$orig.ident))
HNSCC$Sample <- as.factor(as.character(HNSCC$Sample))

###Ranking.
Idents(HNSCC) <- HNSCC$Celltype
HNSCC <- RenameIdents(HNSCC, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "pDCs" = "pDCs", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Epithelial cells" = "Epithelial cells")
table(Idents(HNSCC))
HNSCC$Celltype <- Idents(HNSCC)
table(HNSCC$Celltype)
saveRDS(HNSCC, file = "HNSCC_resolution0.3_annotation_sub.rds")

###Visualization.
#UMAP plot for Figure 1C.
pdf(file="HNSCC_umap_Celltype.pdf",width=5.8, height=4)
DimPlot(HNSCC, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Figure 1D.
Idents(HNSCC) <- 'Celltype'
My_levels <- c( "T cells", "B cells", "Plasma cells", "pDCs", "Mast cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Epithelial cells")
Idents(HNSCC) <- factor(Idents(HNSCC), levels= My_levels)
features = c("KRT19", "EPCAM", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "TCL1A", "LILRA4", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="HNSCC_features.pdf",width=8, height=9)
DotPlot(HNSCC, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2A.
pdf(file="HNSCC_umap_orig.ident.pdf",width=5.6, height=4)
DimPlot(HNSCC, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Figure 1E, F.
###Cluster the Fibroblasts in Primary tumor and Metastatic tumor.
HNSCC_Fib <- subset(HNSCC, subset = Celltype %in% c("Fibroblasts"))
HNSCC_Fib$seurat_clusters <- as.factor(as.character(HNSCC_Fib$seurat_clusters))
HNSCC_Fib$Celltype <- as.factor(as.character(HNSCC_Fib$Celltype))
HNSCC_Fib$orig.ident <- as.factor(as.character(HNSCC_Fib$orig.ident))
saveRDS(HNSCC_Fib, file = "HNSCC_Fib.rds")

#CAFs.
HNSCC_Fib_T <- subset(HNSCC_Fib, subset = Sample %in% c("Primary", "Metastatic"))
HNSCC_Fib_T$seurat_clusters <- as.factor(as.character(HNSCC_Fib_T$seurat_clusters))
HNSCC_Fib_T$Sample <- as.factor(as.character(HNSCC_Fib_T$Sample))
HNSCC_Fib_T$orig.ident <- as.factor(as.character(HNSCC_Fib_T$orig.ident))
Fib <- HNSCC_Fib_T
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
saveRDS(Fib, file = "HNSCC_Fib_T_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4A.
library(clustree)
library(patchwork)
pdf(file = "clustree_HNSCC_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
Fib <- FindClusters(Fib, resolution = 0.2)
saveRDS(Fib, file = "HNSCC_Fib_T_res0.2.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="HNSCC_Fib_T_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="HNSCC_Fib_T_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="HNSCC_Fib_T_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="HNSCC_Fib_T_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5")
pdf(file="HNSCC_Fib_T_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(4)
iCAFs = c(2)
myCAFs = c(3)
eCAFs = c(0, 1, 5)
unknown = c(6)
current.cluster.ids <- c(apCAFs,
                         iCAFs,
                         myCAFs,
                         eCAFs,
                         unknown)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
                     rep("iCAFs",length(iCAFs)),
                     rep("myCAFs",length(myCAFs)),
                     rep("eCAFs",length(eCAFs)),
                     rep("unknown",length(unknown))
)

Fib@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Fib@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Fib$celltype)
table(Fib$seurat_clusters)
saveRDS(Fib, file = "HNSCC_Fib_T_res0.2_annotation.rds")

###Exclude unknown cell subpopulations.
Fib <- subset(Fib, subset = celltype %in% c("apCAFs", "iCAFs", "myCAFs", "eCAFs"))
Fib$seurat_clusters <- as.factor(as.character(Fib$seurat_clusters))
Fib$celltype <- as.factor(as.character(Fib$celltype))
Fib$orig.ident <- as.factor(as.character(Fib$orig.ident))
Fib$Sample <- as.factor(as.character(Fib$Sample))

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "HNSCC_Fib_T_res0.2_annotation_sub.rds")

###Visualization
#UMAP plot for Figure 1E.
pdf(file="HNSCC_Fib_T_umap_celltype.pdf",width=5.3, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Figure 1F.
Idents(Fib) <- Fib$celltype
F_features = c("CXCL14", "CXCL12", "CCL2", "CD34", "IL6", "POSTN", "LOXL2", "MMP14", "RGS5", "ACTA2", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "CD74")
pdf(file="HNSCC_Fib_T_F_features_res0.1_annotation.pdf",width=8, height=6.5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
