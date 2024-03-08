###codes for Supplementary Figure 5C-D.

#1 Supplementary Figure 5C.
#CC.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
folders=list.files('./','^GSM')
folders
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList

###Create one Seurat object
Cervical <- merge(sceList[[1]], 
                  y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                        sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]]), 
                  add.cell.ids = folders, 
                  project = "Cervical_Cancer")
Cervical
table(Cervical$orig.ident)
head(colnames(Cervical))
head(Cervical@meta.data)

###Annotate the Sample.
Idents(Cervical) <- Cervical$orig.ident
Cervical <- RenameIdents(Cervical, "GSM6360680" = "Normal", "GSM6360681" = "Normal", "GSM6360682" = "Normal", "GSM6360683" = "Normal", "GSM6360684" = "HSIL", "GSM6360685" = "HSIL", "GSM6360686" = "Tumor", "GSM6360687" = "Tumor", "GSM6360688" = "Tumor")
table(Idents(Cervical))
Cervical$Sample <- Idents(Cervical)
table(Cervical$Sample)

###Rename the orig.ident.
Idents(Cervical) <- Cervical$orig.ident
Cervical <- RenameIdents(Cervical, "GSM6360680" = "N1_NO_HPV", "GSM6360681" = "N2_NO_HPV", "GSM6360682" = "N1_HPV", "GSM6360683" = "N2_HPV", "GSM6360684" = "HSIL1_HPV", "GSM6360685" = "HSIL2_HPV", "GSM6360686" = "CA1_HPV", "GSM6360687" = "CA2_HPV", "GSM6360688" = "CA3_HPV")
table(Idents(Cervical))
Cervical$orig.ident <- Idents(Cervical)
table(Cervical$orig.ident)
head(Cervical@meta.data)
saveRDS(Cervical, file = "Cervical_Seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(Cervical, expression = nFeature_RNA > 200)
selected_f <- rownames(Cervical)[Matrix::rowSums(Cervical) > 3]
Cervical <- subset(Cervical, features = selected_f, cells = selected_c)
dim(Cervical)

###QC metrics and filtering.
Cervical[["percent.mt"]] <- PercentageFeatureSet(Cervical, pattern = "^MT-")
Idents(Cervical) <- Cervical$orig.ident
pdf(file="Cervical.featureViolin.pdf",width=9,height=6)           
VlnPlot(object = Cervical, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="Cervical.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(Cervical, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cervical, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
Cervical <- subset(Cervical, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(Cervical)

###Visualization.
#VlnPlot for Supplementary Figure 1D.
Idents(Cervical) <- Cervical$orig.ident
pdf(file="Cervical.featureViolin_QC.pdf",width=8,height=6)           
VlnPlot(object = Cervical, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
Cervical <- NormalizeData(Cervical, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
Cervical <- FindVariableFeatures(Cervical, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(Cervical)
Cervical <- ScaleData(Cervical, features = all.genes)
Cervical <- ScaleData(Cervical, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
Cervical <- RunPCA(Cervical, features = VariableFeatures(object = Cervical))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
Cervical$orig.ident <- as.factor(as.character(Cervical$orig.ident))
Cervical <- RunHarmony(Cervical, group.by.vars = "orig.ident")
Cervical <- RunUMAP(Cervical, reduction = "harmony", dims = 1:20)
Cervical <- RunTSNE(Cervical, reduction = "harmony", dims = 1:20)
names(Cervical@reductions)

###Cluster the cells.
Cervical <- FindNeighbors(Cervical, reduction = "harmony", dims = 1:20)
saveRDS(Cervical, file = "Cervical_before_resolution.rds")

#Determine resolution.
test <- Cervical
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3D.
library(clustree)
library(patchwork)
pdf(file = "clustree_Cervical.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.3.
Cervical <- FindClusters(Cervical, resolution = 0.3)
saveRDS(Cervical, file = "Cervical_resolution_0.3.rds")

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
pdf(file="Cervical_marker.pdf",width=44, height=33)
VlnPlot(Cervical, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="Cervical_NK.pdf",width=44, height=14)
VlnPlot(Cervical, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="Cervical_Epithelial.pdf",width=44, height=14)
VlnPlot(Cervical, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="Cervical_Plasma_cells.pdf",width=44, height=14)
VlnPlot(Cervical, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="Cervical_B_marker.pdf",width=44, height=19)
VlnPlot(Cervical, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="Cervical_Fibroblasts.pdf",width=44, height=14)
VlnPlot(Cervical, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="Cervical_Endothelial.pdf",width=44, height=5)
VlnPlot(Cervical, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="Cervical_T_marker.pdf",width=44, height=14)
VlnPlot(Cervical, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="Cervical_Mast.pdf",width=44, height=5)
VlnPlot(Cervical, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="Cervical_myeloid.pdf",width=44, height=9.5)
VlnPlot(Cervical, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="Cervical_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(Cervical, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="Cervical_pDC.pdf",width=44, height=20)
VlnPlot(Cervical, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="Cervical_neutrophil.pdf",width=44, height=5)
VlnPlot(Cervical, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="Cervical_melanocyte.pdf",width=44, height=30)
VlnPlot(Cervical, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="Cervical_clear_cell.pdf",width=44, height=10)
VlnPlot(Cervical, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(6, 14)
Myeloid_cells = c(2, 3)
Mast_cells = c(11)
NK_cells = c(8)
B_cells = c(10)
Plasma_cells = c(5)
T_cells = c(0, 1, 13)
Epithelial_cells = c(4, 7, 9, 15)
Endothelial_cells = c(12)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Mast_cells,
                         NK_cells,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Epithelial_cells,
                         Endothelial_cells)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("NK_cells",length(NK_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells))
)

Cervical@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(Cervical@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Cervical$Celltype)
table(Cervical$seurat_clusters)

###Rename.
Idents(Cervical) <- Cervical$Celltype
Cervical <- RenameIdents(Cervical, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells", "Plasma_cells" = "Plasma cells", "NK_cells" = "NK cells")
table(Idents(Cervical))
Cervical$Celltype <- Idents(Cervical)
table(Cervical$Celltype)
saveRDS(Cervical, file = "Cervical_resolution0.3_annotation.rds")

###Ranking.
Idents(Cervical) <- Cervical$Celltype
Cervical <- RenameIdents(Cervical, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "NK cells" = "NK cells", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Epithelial cells" = "Epithelial cells")
table(Idents(Cervical))
Cervical$Celltype <- Idents(Cervical)
table(Cervical$Celltype)
saveRDS(Cervical, file = "Cervical_resolution0.3_annotation_2.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5C: left.
pdf(file="Cervical_umap_Celltype.pdf",width=5.8, height=4)
DimPlot(Cervical, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5C: right.
Idents(Cervical) <- 'Celltype'
My_levels <- c( "T cells", "B cells", "Plasma cells", "NK cells", "Mast cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Epithelial cells")
Idents(Cervical) <- factor(Idents(Cervical), levels= My_levels)
features = c("KRT19", "EPCAM", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "GNLY", "KLRD1", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="Cervical_features.pdf",width=8, height=9.5)
DotPlot(Cervical, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2D.
pdf(file="Cervical_umap_orig.ident.pdf",width=5.6, height=4)
DimPlot(Cervical, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5D.
###Cluster the Fibroblasts in tumor.
Cervical_Fib <- subset(Cervical, subset = Celltype %in% c("Fibroblasts"))
Cervical_Fib$seurat_clusters <- as.factor(as.character(Cervical_Fib$seurat_clusters))
Cervical_Fib$Celltype <- as.factor(as.character(Cervical_Fib$Celltype))
Cervical_Fib$orig.ident <- as.factor(as.character(Cervical_Fib$orig.ident))
saveRDS(Cervical_Fib, file = "Cervical_Fib.rds")

#CAFs.
Cervical_Fib_T <- subset(Cervical_Fib, subset = Sample %in% c("Tumor"))
Cervical_Fib_T$seurat_clusters <- as.factor(as.character(Cervical_Fib_T$seurat_clusters))
Cervical_Fib_T$Sample <- as.factor(as.character(Cervical_Fib_T$Sample))
Cervical_Fib_T$orig.ident <- as.factor(as.character(Cervical_Fib_T$orig.ident))
Fib <- Cervical_Fib_T
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
saveRDS(Fib, file = "Cervical_Fib_T_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4D.
library(clustree)
library(patchwork)
pdf(file = "clustree_Cervical_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 1.
Fib <- FindClusters(Fib, resolution = 1)
saveRDS(Fib, file = "Cervical_Fib_T_res1.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="Cervical_Fib_T_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="Cervical_Fib_T_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="Cervical_Fib_T_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="Cervical_Fib_T_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5")
pdf(file="Cervical_Fib_T_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(2)
iCAFs = c(0)
myCAFs = c(4)
eCAFs = c(5)
unknown_CAFs = c(1, 3)
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
saveRDS(Fib, file = "Cervical_Fib_T_res1_annotation.rds")

###Exclude unknown cell subpopulations.
Fib <- subset(Fib, subset = celltype %in% c("apCAFs", "myCAFs", "eCAFs", "iCAFs"))
Fib$seurat_clusters <- as.factor(as.character(Fib$seurat_clusters))
Fib$celltype <- as.factor(as.character(Fib$celltype))
Fib$orig.ident <- as.factor(as.character(Fib$orig.ident))

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "Cervical_Fib_T_res1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5D: left.
pdf(file="Cervical_Fib_T_umap_celltype.pdf",width=5.5, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5D: right.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "eCAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "POSTN", "MMP14", "RGS5", "ACTA2", "HLA-DPB1", "HLA-DRA", "CD74")
pdf(file="Cervical_Fib_T_F_features_res1_annotation.pdf",width=8, height=5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
