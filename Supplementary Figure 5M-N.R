###codes for Supplementary Figure 5M-N.

###Export the raw count matrix in scanpy
#linux
#conda activate scanpy
#python3
#import numpy as np
#import pandas as pd
#import scanpy as sc

#sc.settings.verbosity = 3             # 设置日志等级: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_header()
#sc.settings.set_figure_params(dpi=80, facecolor='white')

#results_file = '/media/Cambyses/Entertainment/RCC/RCC.h5ad'
#RCC=sc.read('/media/Cambyses/Entertainment/RCC/RCC_upload_final_raw_counts.h5ad')
#RCC
#import scvelo as scv
#scv.DataFrame(RCC.X)
#RCC.obs
#RCC.var
#pd.DataFrame(data=RCC.X.todense().T, index=RCC.var_names,columns=RCC.obs_names).to_csv('RCC_matrix.csv',sep="\t",float_format='%.0f') 
#gzip RCC_matrix.csv

#R
#1 Supplementary Figure 5M.
#RCC.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
library(Seurat)
library(scuttle)
library(ggplot2)
library(scales)
library(ggsci)
t <- readSparseCounts("RCC_matrix.csv.gz")
rcc <- CreateSeuratObject(counts = t, project = "RCC")
rcc
head(rcc@meta.data)
table(rcc$orig.ident)

#Keep only samples of Normal kidney, Primary tumor, and Tumor-normal interface.
RCC <- subset(rcc, subset=orig.ident %in% c("5739STDY7958793", "5739STDY7958794", "5739STDY7958795", "5739STDY7958805", "5739STDY7958806", "5739STDY7958807", "5739STDY7958808", "5739STDY7958809", "5739STDY7958811", "5739STDY7958812", "5739STDY8351215", "5739STDY8351216", "5739STDY8351217", "5739STDY8351218", "5739STDY8351221", "5739STDY8351222", "5739STDY8351223", "5739STDY8351224", "5739STDY8351225", "5739STDY8351226", "5739STDY8351229", "5739STDY8351232", "5739STDY8351233", "5739STDY8351234", "5739STDY8351240", "5739STDY8351241", "5739STDY8351242", "5739STDY8351249", "5739STDY8351250", "5739STDY8351256", "5739STDY8351257", "5739STDY8351264", "5739STDY8351266", "5739STDY9266968", "5739STDY9266969", "5739STDY9266970", "5739STDY9266971", "5739STDY9266973", "5739STDY9266974", "5739STDY9266976", "5739STDY9266977", "5739STDY9266978", "5739STDY9266979", "5739STDY9266980", "5739STDY9266981", "5739STDY9266985", "5739STDY9266990", "5739STDY9266992", "5739STDY9266993", "5739STDY9266994", "5739STDY9266995", "5739STDY9266997", "5739STDY9266998"))
RCC$orig.ident <- as.factor(as.character(RCC$orig.ident))
RCC$orig.ident <- paste0("X", RCC$orig.ident)

#Rename cells name.
head(x = colnames(x = RCC))
RCC <- RenameCells(object = RCC, new.names = paste0("X", colnames(x = RCC)))
head(x = colnames(x = RCC))
name <- gsub("-", ".", colnames(x = RCC))
RCC <- RenameCells(object = RCC, new.names = name)
head(x = colnames(x = RCC))
head(RCC@meta.data)
table(RCC$orig.ident)

#Comment sample information.
library(dplyr)
RCC@meta.data <- RCC@meta.data %>%
  mutate(
    Sample = case_when(
      orig.ident %in% c("X5739STDY7958793", "X5739STDY7958805", "X5739STDY7958806", "X5739STDY7958807", "X5739STDY7958808", "X5739STDY7958809", "X5739STDY8351215", "X5739STDY8351216", "X5739STDY8351217", "X5739STDY8351218", "X5739STDY8351221", "X5739STDY8351222", "X5739STDY8351223", "X5739STDY8351224", "X5739STDY8351225", "X5739STDY8351226", "X5739STDY8351232", "X5739STDY8351233", "X5739STDY8351234", "X5739STDY8351240", "X5739STDY8351241", "X5739STDY8351242", "X5739STDY9266968", "X5739STDY9266969", "X5739STDY9266970", "X5739STDY9266971", "X5739STDY9266976", "X5739STDY9266977", "X5739STDY9266978", "X5739STDY9266979", "X5739STDY9266985", "X5739STDY9266990", "X5739STDY9266992", "X5739STDY9266993", "X5739STDY9266994", "X5739STDY9266995") ~ "Primary tumor",
      orig.ident %in% c("X5739STDY7958795", "X5739STDY7958811", "X5739STDY8351249", "X5739STDY8351250", "X5739STDY8351256", "X5739STDY9266973", "X5739STDY9266980", "X5739STDY9266997") ~ "Tumor-normal interface",
      TRUE ~ "Normal kidney"  #If all conditions are not met, set it to "Normal kidney"
    )
  )
table(RCC$Sample)
table(RCC$orig.ident)
head(RCC@meta.data)
saveRDS(RCC, file = "RCC_Seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(RCC, expression = nFeature_RNA > 200)
selected_f <- rownames(RCC)[Matrix::rowSums(RCC) > 3]
RCC <- subset(RCC, features = selected_f, cells = selected_c)
dim(RCC)

###QC metrics and filtering.
RCC[["percent.mt"]] <- PercentageFeatureSet(RCC, pattern = "^MT-")
Idents(RCC) <- RCC$orig.ident
pdf(file="RCC.featureViolin.pdf",width=36,height=6)           
VlnPlot(object = RCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="RCC.FeatureScatter.pdf",width=24,height=6)
plot1 <- FeatureScatter(RCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
RCC <- subset(RCC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(RCC)

###Visualization.
#VlnPlot for Supplementary Figure 1I.
Idents(RCC) <- RCC$orig.ident
pdf(file="RCC.featureViolin_QC.pdf",width=35,height=6)           
VlnPlot(object = RCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
RCC <- NormalizeData(RCC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
RCC <- FindVariableFeatures(RCC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(RCC)
RCC <- ScaleData(RCC, features = all.genes)
RCC <- ScaleData(RCC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
RCC <- RunPCA(RCC, features = VariableFeatures(object = RCC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
RCC <- RunHarmony(RCC, group.by.vars = "orig.ident")
RCC <- RunUMAP(RCC, reduction = "harmony", dims = 1:20)
RCC <- RunTSNE(RCC, reduction = "harmony", dims = 1:20)
names(RCC@reductions)

###Cluster the cells.
RCC <- FindNeighbors(RCC, reduction = "harmony", dims = 1:20)
saveRDS(RCC, file = "RCC_before_resolution.rds")

#Determine resolution.
test <- RCC
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3I.
library(clustree)
library(patchwork)
pdf(file = "clustree_RCC.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
RCC <- FindClusters(RCC, resolution = 0.1)
saveRDS(RCC, file = "RCC_resolution_0.1.rds")

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
pdf(file="RCC_marker.pdf",width=44, height=33)
VlnPlot(RCC, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="RCC_NK.pdf",width=44, height=14)
VlnPlot(RCC, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="RCC_Epithelial.pdf",width=44, height=14)
VlnPlot(RCC, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="RCC_Plasma_cells.pdf",width=44, height=14)
VlnPlot(RCC, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="RCC_B_marker.pdf",width=44, height=19)
VlnPlot(RCC, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="RCC_Fibroblasts.pdf",width=44, height=14)
VlnPlot(RCC, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="RCC_Endothelial.pdf",width=44, height=5)
VlnPlot(RCC, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="RCC_T_marker.pdf",width=44, height=14)
VlnPlot(RCC, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="RCC_Mast.pdf",width=44, height=5)
VlnPlot(RCC, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="RCC_myeloid.pdf",width=44, height=9.5)
VlnPlot(RCC, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="RCC_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(RCC, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="RCC_pDC.pdf",width=44, height=20)
VlnPlot(RCC, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="RCC_neutrophil.pdf",width=44, height=5)
VlnPlot(RCC, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="RCC_melanocyte.pdf",width=44, height=30)
VlnPlot(RCC, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="RCC_clear_cell.pdf",width=44, height=10)
VlnPlot(RCC, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(8)
Myeloid_cells = c(3, 9)
Mast_cells = c(10)
B_cells = c(7)
NK_cells = c(2)
T_cells = c(0, 1, 5)
Clear_cells = c(4)
Endothelial_cells = c(6)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Mast_cells,
                         B_cells,
                         NK_cells,
                         T_cells,
                         Clear_cells,
                         Endothelial_cells)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Clear_cells",length(Clear_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells))
)

RCC@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(RCC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(RCC$Celltype)

###Rename.
Idents(RCC) <- RCC$Celltype
RCC <- RenameIdents(RCC, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "B_cells" = "B cells", "Clear_cells" = "Clear cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells", "NK_cells" = "NK cells")
table(Idents(RCC))
RCC$Celltype <- Idents(RCC)
table(RCC$Celltype)
saveRDS(RCC, file = "RCC_resolution0.1_annotation.rds")

###Ranking.
Idents(RCC) <- RCC$Celltype
RCC <- RenameIdents(RCC, "T cells" = "T cells", "B cells" = "B cells", "NK cells" = "NK cells", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Clear cells" = "Clear cells")
table(Idents(RCC))
RCC$Celltype <- Idents(RCC)
table(RCC$Celltype)
saveRDS(RCC, file = "RCC_resolution0.1_annotation2.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5M: left.
pdf(file="RCC_umap_Celltype.pdf",width=5.6, height=4)
DimPlot(RCC, reduction = "umap",label=T, group.by = "Celltype", repel = T, raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5M: right.
Idents(RCC) <- 'Celltype'
My_levels <- c( "T cells", "B cells", "NK cells", "Mast cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Clear cells")
Idents(RCC) <- factor(Idents(RCC), levels= My_levels)
features = c("GATM", "CA9", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "FGFBP2", "KLRD1", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="RCC_features.pdf",width=8, height=8.5)
DotPlot(RCC, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2I.
pdf(file="RCC_umap_orig.ident.pdf",width=10, height=4)
DimPlot(RCC, reduction = "umap",label=T, group.by = "orig.ident", raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5N.
###Cluster the Fibroblasts in Primary tumor.
RCC_Fib <- subset(RCC, subset = Celltype %in% c("Fibroblasts"))
RCC_Fib$seurat_clusters <- as.factor(as.character(RCC_Fib$seurat_clusters))
RCC_Fib$Celltype <- as.factor(as.character(RCC_Fib$Celltype))
RCC_Fib$orig.ident <- as.factor(as.character(RCC_Fib$orig.ident))
saveRDS(RCC_Fib, file = "RCC_Fib.rds")

#CAFs.
RCC_Fib_T <- subset(RCC_Fib, subset = Sample %in% c("Primary tumor"))
RCC_Fib_T$seurat_clusters <- as.factor(as.character(RCC_Fib_T$seurat_clusters))
RCC_Fib_T$Celltype <- as.factor(as.character(RCC_Fib_T$Celltype))
RCC_Fib_T$orig.ident <- as.factor(as.character(RCC_Fib_T$orig.ident))
Fib <- RCC_Fib_T
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
saveRDS(Fib, file = "RCC_Fib_T_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4I.
library(clustree)
library(patchwork)
pdf(file = "clustree_RCC_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.3.
Fib <- FindClusters(Fib, resolution = 0.3)
saveRDS(Fib, file = "RCC_Fib_T_res0.3.rds")

###classic celltype marker in clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="RCC_Fib_T_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="RCC_Fib_T_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="RCC_Fib_T_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="RCC_Fib_T_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5")
pdf(file="RCC_Fib_T_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(6)  #There are two CAFs clusters with high expression of CD74 and HLA-D genes, and we only selected one CAFs cluster 6 with higher expression of CD74 as apCAFs.
myCAFs = c(0, 4)
eCAFs = c(7)
unknown_CAFs = c(1, 2, 3, 5)
current.cluster.ids <- c(apCAFs,
                         myCAFs,
                         eCAFs,
                         unknown_CAFs)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
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
saveRDS(Fib, file = "RCC_Fib_T_res0.3_annotation.rds")

###Exclude unknown cell subpopulations.
Fib <- subset(Fib, subset = celltype %in% c("apCAFs", "myCAFs", "eCAFs"))
Fib$seurat_clusters <- as.factor(as.character(Fib$seurat_clusters))
Fib$celltype <- as.factor(as.character(Fib$celltype))
Fib$orig.ident <- as.factor(as.character(Fib$orig.ident))

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "RCC_Fib_T_res0.3_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5N: left.
pdf(file="RCC_Fib_T_umap_celltype.pdf",width=5.3, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5N: right.
Idents(Fib) <- Fib$celltype
F_features = c("POSTN", "MMP14", "RGS5", "ACTA2", "HLA-DPB1", "HLA-DPA1", "CD74")
pdf(file="RCC_Fib_T_F_features_res0.3_annotation.pdf",width=8, height=4.5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
