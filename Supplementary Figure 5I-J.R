###codes for Supplementary Figure 5I-J.

#1 Supplementary Figure 5I.
#AM.
###Read and preprocess the original count matrix.
##1
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708993")
folders=list.files('./','^SRR')
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
GSM5708993 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708993
table(GSM5708993$orig.ident)
head(colnames(GSM5708993))
##2
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708994")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708994 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708994
table(GSM5708994$orig.ident)
head(colnames(GSM5708994))

##3
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708995")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708995 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708995
table(GSM5708995$orig.ident)
head(colnames(GSM5708995))

##4
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708996")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708996 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708996
table(GSM5708996$orig.ident)
head(colnames(GSM5708996))

##5
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708997")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708997 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708997
table(GSM5708997$orig.ident)
head(colnames(GSM5708997))

##6
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708998")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708998 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708998
table(GSM5708998$orig.ident)
head(colnames(GSM5708998))

##7
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5708999")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5708999 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5708999
table(GSM5708999$orig.ident)
head(colnames(GSM5708999))

##8
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5709000")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5709000 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5709000
table(GSM5709000$orig.ident)
head(colnames(GSM5709000))

##9
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889/GSM5709001")
folders=list.files('./','^SRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
GSM5709001 <- merge(sceList[[1]], 
                    y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                          sceList[[6]],sceList[[7]],sceList[[8]]), 
                    add.cell.ids = folders, 
                    project = "AM")
GSM5709001
table(GSM5709001$orig.ident)
head(colnames(GSM5709001))

###annotation
GSM5708993$Sample = "metastatic"
GSM5708994$Sample = "primary"
GSM5708995$Sample = "primary"
GSM5708996$Sample = "primary"
GSM5708997$Sample = "metastatic"
GSM5708998$Sample = "primary"
GSM5708999$Sample = "metastatic"
GSM5709000$Sample = "metastatic"
GSM5709001$Sample = "primary"

###Rename the orig.ident.
GSM5708993$orig.ident = "AM1"
GSM5708994$orig.ident = "AM2"
GSM5708995$orig.ident = "AM3"
GSM5708996$orig.ident = "AM4"
GSM5708997$orig.ident = "AM5"
GSM5708998$orig.ident = "AM6"
GSM5708999$orig.ident = "AM7"
GSM5709000$orig.ident = "AM8M"
GSM5709001$orig.ident = "AM8P"

###Create one Seurat object
AM  <- merge(GSM5708993, y = c(GSM5708994, GSM5708995, GSM5708996, GSM5708997, GSM5708998, GSM5708999, GSM5709000, GSM5709001), add.cell.ids = c("AM93", "AM94", "AM95", "AM96", "AM97", "AM98", "AM99", "AM00", "AM01"), project = "AM")
AM
head(colnames(AM))
tail(colnames(AM))
table(AM$orig.ident)
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
saveRDS(AM, file = "AM_seurat.rds")

###Preliminary filtering.
selected_c <- WhichCells(AM, expression = nFeature_RNA > 200)
selected_f <- rownames(AM)[Matrix::rowSums(AM) > 3]
AM <- subset(AM, features = selected_f, cells = selected_c)
dim(AM)

###QC metrics and filtering.
AM[["percent.mt"]] <- PercentageFeatureSet(AM, pattern = "^MT-")
Idents(AM) <- AM$orig.ident
pdf(file="AM.featureViolin.pdf",width=12,height=6)           
VlnPlot(object = AM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="AM.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(AM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
AM <- subset(AM, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(AM)

###Visualization.
#VlnPlot for Supplementary Figure 1H.
Idents(AM) <- AM$orig.ident
pdf(file="AM.featureViolin_QC.pdf",width=8,height=6)           
VlnPlot(object = AM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
AM <- NormalizeData(AM, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
AM <- FindVariableFeatures(AM, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(AM)
AM <- ScaleData(AM, features = all.genes)
AM <- ScaleData(AM, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
AM <- RunPCA(AM, features = VariableFeatures(object = AM))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
AM <- RunHarmony(AM, group.by.vars = "orig.ident")
AM <- RunUMAP(AM, reduction = "harmony", dims = 1:20)
AM <- RunTSNE(AM, reduction = "harmony", dims = 1:20)
names(AM@reductions)

###Cluster the cells.
AM <- FindNeighbors(AM, reduction = "harmony", dims = 1:20)
saveRDS(AM, file = "AM_before_resolution.rds")

#Determine resolution.
test <- AM
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3G.
library(clustree)
library(patchwork)
pdf(file = "clustree_AM.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
AM <- FindClusters(AM, resolution = 0.1)
saveRDS(AM, file = "AM_resolution_0.1.rds")

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
pdf(file="AM_marker.pdf",width=44, height=33)
VlnPlot(AM, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="AM_NK.pdf",width=44, height=14)
VlnPlot(AM, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="AM_Epithelial.pdf",width=44, height=14)
VlnPlot(AM, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="AM_Plasma_cells.pdf",width=44, height=14)
VlnPlot(AM, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="AM_B_marker.pdf",width=44, height=19)
VlnPlot(AM, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="AM_Fibroblasts.pdf",width=44, height=14)
VlnPlot(AM, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="AM_Endothelial.pdf",width=44, height=5)
VlnPlot(AM, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="AM_T_marker.pdf",width=44, height=14)
VlnPlot(AM, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="AM_Mast.pdf",width=44, height=5)
VlnPlot(AM, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="AM_myeloid.pdf",width=44, height=9.5)
VlnPlot(AM, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="AM_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(AM, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="AM_pDC.pdf",width=44, height=20)
VlnPlot(AM, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="AM_neutrophil.pdf",width=44, height=5)
VlnPlot(AM, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="AM_melanocyte.pdf",width=44, height=30)
VlnPlot(AM, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="AM_clear_cell.pdf",width=44, height=10)
VlnPlot(AM, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(11)
Myeloid_cells = c(5)
Plasma_cells = c(13)
T_cells = c(1, 14)
Melanocytes = c(0, 2, 3, 4, 6, 7, 10, 12, 16)
Endothelial_cells = c(8, 15)
pDCs = c(17)
unknown = c(9)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Plasma_cells,
                         T_cells,
                         Melanocytes,
                         Endothelial_cells,
                         pDCs,
                         unknown)

new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Melanocytes",length(Melanocytes)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("pDCs",length(pDCs)),
                     rep("unknown",length(unknown))
)

AM@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(AM@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(AM@meta.data$Celltype)
table(AM$seurat_clusters)

###Rename.
Idents(AM) <- AM$Celltype
AM <- RenameIdents(AM, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "Plasma_cells" = "Plasma cells", "Myeloid_cells" = "Myeloid cells")
table(Idents(AM))
AM$Celltype <- Idents(AM)
table(AM$Celltype)
saveRDS(AM, file = "AM_resolution0.1_annotation.rds")

###Exclude unknown cell subpopulations.
AM <- subset(AM, subset = Celltype %in% c("Fibroblasts", "Myeloid cells", "pDCs", "Plasma cells", "T cells", "Melanocytes", "Endothelial cells"))
AM$seurat_clusters <- as.factor(as.character(AM$seurat_clusters))
AM$Celltype <- as.factor(as.character(AM$Celltype))
AM$orig.ident <- as.factor(as.character(AM$orig.ident))
AM$Sample <- as.factor(as.character(AM$Sample))

###Ranking.
Idents(AM) <- AM$Celltype
AM <- RenameIdents(AM, "T cells" = "T cells", "Plasma cells" = "Plasma cells", "pDCs" = "pDCs", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Melanocytes" = "Melanocytes")
table(Idents(AM))
AM$Celltype <- Idents(AM)
table(AM$Celltype)
saveRDS(AM, file = "AM_resolution0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5I: left.
Idents(AM) <- AM$Celltype
pdf(file="AM_umap_Celltype.pdf",width=5.7, height=4)
DimPlot(AM, reduction = "umap",label=T, group.by = "Celltype", repel = T, raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5I: right.
Idents(AM) <- 'Celltype'
My_levels <- c("T cells", "Plasma cells", "pDCs", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Melanocytes")
Idents(AM) <- factor(Idents(AM), levels= My_levels)
features = c("MITF", "MLANA", "PMEL", "TYR", "VWF", "PECAM1", "COL3A1", "COL1A1", "AIF1", "LYZ", "PTGDS", "LILRA4", "IGHG1", "JCHAIN", "CD3E", "CD3D")
pdf(file="AM_features.pdf",width=8, height=8.5)
DotPlot(AM, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2G.
pdf(file="AM_umap_orig.ident.pdf",width=5, height=4)
DimPlot(AM, reduction = "umap",label=T, group.by = "orig.ident", raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5J.
###Cluster the Fibroblasts in tumor.
AM_Fib <- subset(AM, subset = Celltype %in% c("Fibroblasts"))
AM_Fib$Celltype <- as.factor(as.character(AM_Fib$Celltype))
AM_Fib$orig.ident <- as.factor(as.character(AM_Fib$orig.ident))
saveRDS(AM_Fib, file = "AM_Fib.rds")
Fib <- AM_Fib
table(Fib$Celltype)
table(Fib$Sample)

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
saveRDS(Fib, file = "AM_Fib_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4G.
library(clustree)
library(patchwork)
pdf(file = "clustree_AM_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
Fib <- FindClusters(Fib, resolution = 0.1)
saveRDS(Fib, file = "AM_Fib_res0.1.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="AM_Fib_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="AM_Fib_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="AM_Fib_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="AM_Fib_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5", "MYH11")
pdf(file="AM_Fib_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(7)
eCAFs = c(3, 4)
myCAFs = c(1)
MYH11_CAFs = c(2)
iCAFs = c(0, 5, 6)
current.cluster.ids <- c(apCAFs,
                         eCAFs,
                         myCAFs,
                         MYH11_CAFs,
                         iCAFs)

new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
                     rep("eCAFs",length(eCAFs)),
                     rep("myCAFs",length(myCAFs)),
                     rep("MYH11_CAFs",length(MYH11_CAFs)),
                     rep("iCAFs",length(iCAFs))
)

Fib@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Fib@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Fib$celltype)
table(Fib$seurat_clusters)

###Rename.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "MYH11_CAFs" = "MYH11+ CAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "AM_Fib_res0.1_annotation.rds")

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "MYH11+ CAFs" = "MYH11+ CAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "AM_Fib_res0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5J: left.
pdf(file="AM_Fib_umap_celltype.pdf",width=5.6, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5J: right.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "MYH11+ CAFs", "eCAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "CXCL12", "CCL2", "IL6", "POSTN", "LOXL2", "MMP14", "MYH11", "RGS5", "ACTA2", "HLA-DRB1", "HLA-DRA", "CD74")
pdf(file="AM_Fib_F_features_res0.1_annotation.pdf",width=8, height=7)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
