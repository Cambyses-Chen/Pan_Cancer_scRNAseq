###codes for Supplementary Figure 5E-F.

#1 Supplementary Figure 5E.
#CRC.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
folders=list.files('./','^HRR')
folders
library(Seurat)
library(ggplot2)
library(ggsci)
library(scales)
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})
sceList
###Create one Seurat object
CRC <- merge(sceList[[1]], 
             y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                   sceList[[6]],sceList[[7]],sceList[[8]], sceList[[9]], sceList[[10]]), 
             add.cell.ids = folders, 
             project = "CRC")
CRC
head(CRC@meta.data)

###Annotate the Sample.
table(CRC$orig.ident)
CRC$Sample <- ifelse(test = CRC$orig.ident %in% c("HRR273223", "HRR273225", "HRR273227", "HRR273229", "HRR273231"), yes = "Normal", no = "Tumor")
table(CRC$Sample)

###Rename the orig.ident.
Idents(CRC) <- CRC$orig.ident
CRC <- RenameIdents(CRC, "HRR273223" = "P1-CRC-N", "HRR273224" = "P1-CRC-T", "HRR273225" = "P2-CRC-N", "HRR273226" = "P2-CRC-T", "HRR273227" = "P3-CRC-N", "HRR273228" = "P3-CRC-T", "HRR273229" = "P4-CRC-N", "HRR273230" = "P4-CRC-T", "HRR273231" = "P5-CRC-N", "HRR273232" = "P5-CRC-T")
table(Idents(CRC))
CRC$orig.ident <- Idents(CRC)
table(CRC$orig.ident)
head(CRC@meta.data)
saveRDS(CRC, file = "CRC_HRA000979_Merge.rds")

###Preliminary filtering.
selected_c <- WhichCells(CRC, expression = nFeature_RNA > 200)
selected_f <- rownames(CRC)[Matrix::rowSums(CRC) > 3]
CRC <- subset(CRC, features = selected_f, cells = selected_c)
dim(CRC)

###QC metrics and filtering.
CRC[["percent.mt"]] <- PercentageFeatureSet(CRC, pattern = "^MT-")
Idents(CRC) <- CRC$orig.ident
pdf(file="CRC.featureViolin.pdf",width=9,height=6)           
VlnPlot(object = CRC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="CRC.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(CRC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CRC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
CRC <- subset(CRC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(CRC)

###Visualization.
#VlnPlot for Supplementary Figure 1F.
Idents(CRC) <- CRC$orig.ident
pdf(file="CRC.featureViolin_QC.pdf",width=8,height=6)           
VlnPlot(object = CRC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
CRC <- NormalizeData(CRC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
CRC <- FindVariableFeatures(CRC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(CRC)
CRC <- ScaleData(CRC, features = all.genes)
CRC <- ScaleData(CRC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
CRC <- RunPCA(CRC, features = VariableFeatures(object = CRC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
CRC <- RunHarmony(CRC, group.by.vars = "orig.ident")
CRC <- RunUMAP(CRC, reduction = "harmony", dims = 1:20)
CRC <- RunTSNE(CRC, reduction = "harmony", dims = 1:20)
names(CRC@reductions)

###Cluster the cells.
CRC <- FindNeighbors(CRC, reduction = "harmony", dims = 1:20)
saveRDS(CRC, file = "CRC_HRA000979_before_resolution.rds")

#Determine resolution.
test <- CRC
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3E.
library(clustree)
library(patchwork)
pdf(file = "clustree_CRC.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
CRC <- FindClusters(CRC, resolution = 0.1)
saveRDS(CRC, file = "CRC_HRA000979_resolution_0.1.rds")

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
pdf(file="CRC_marker.pdf",width=44, height=33)
VlnPlot(CRC, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="CRC_NK.pdf",width=44, height=14)
VlnPlot(CRC, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="CRC_Epithelial.pdf",width=44, height=14)
VlnPlot(CRC, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="CRC_Plasma_cells.pdf",width=44, height=14)
VlnPlot(CRC, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="CRC_B_marker.pdf",width=44, height=19)
VlnPlot(CRC, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="CRC_Fibroblasts.pdf",width=44, height=14)
VlnPlot(CRC, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="CRC_Endothelial.pdf",width=44, height=5)
VlnPlot(CRC, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="CRC_T_marker.pdf",width=44, height=14)
VlnPlot(CRC, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="CRC_Mast.pdf",width=44, height=5)
VlnPlot(CRC, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="CRC_myeloid.pdf",width=44, height=9.5)
VlnPlot(CRC, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="CRC_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(CRC, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="CRC_pDC.pdf",width=44, height=20)
VlnPlot(CRC, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="CRC_neutrophil.pdf",width=44, height=5)
VlnPlot(CRC, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="CRC_melanocyte.pdf",width=44, height=30)
VlnPlot(CRC, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="CRC_clear_cell.pdf",width=44, height=10)
VlnPlot(CRC, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(4, 9, 11)
Myeloid_cells = c(7)
B_cells = c(5)
Plasma_cells = c(1)
T_cells = c(0, 2)
Epithelial_cells = c(3)
Endothelial_cells = c(8)
Mast_cells = c(6)
unknown = c(10)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Epithelial_cells,
                         Endothelial_cells,
                         Mast_cells,
                         unknown)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("unknown",length(unknown))
)

CRC@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(CRC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(CRC@meta.data$Celltype)
table(CRC$seurat_clusters)

###Rename.
Idents(CRC) <- CRC$Celltype
CRC <- RenameIdents(CRC, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "Plasma_cells" = "Plasma cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells")
table(Idents(CRC))
CRC$Celltype <- Idents(CRC)
table(CRC$Celltype)
saveRDS(CRC, file = "CRC_HRA000979_resolution0.1_annotation.rds")

###Exclude unknown cell subpopulations.
CRC <- subset(CRC, subset = Celltype %in% c("Fibroblasts", "Myeloid cells", "B cells", "Plasma cells", "T cells", "Epithelial cells", "Endothelial cells", "Mast cells"))
CRC$seurat_clusters <- as.factor(as.character(CRC$seurat_clusters))
CRC$Celltype <- as.factor(as.character(CRC$Celltype))
CRC$orig.ident <- as.factor(as.character(CRC$orig.ident))
CRC$Sample <- as.factor(as.character(CRC$Sample))

###Ranking.
Idents(CRC) <- CRC$Celltype
CRC <- RenameIdents(CRC, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Epithelial cells" = "Epithelial cells")
table(Idents(CRC))
CRC$Celltype <- Idents(CRC)
table(CRC$Celltype)
saveRDS(CRC, file = "CRC_HRA000979_resolution0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5E: left.
pdf(file="CRC_HRA000979_umap_Celltype.pdf",width=5.8, height=4)
DimPlot(CRC, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5E: right.
Idents(CRC) <- 'Celltype'
My_levels <- c("T cells", "B cells", "Plasma cells", "Mast cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Epithelial cells")
Idents(CRC) <- factor(Idents(CRC), levels= My_levels)
features = c("KRT19", "EPCAM", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="CRC_HRA000979_features.pdf",width=8, height=8.5)
DotPlot(CRC, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2E.
pdf(file="CRC_umap_orig.ident.pdf",width=5.5, height=4)
DimPlot(CRC, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5F.
###Cluster the Fibroblasts in tumor.
CRC_Fib <- subset(CRC, subset = Celltype %in% c("Fibroblasts"))
CRC_Fib$Celltype <- as.factor(as.character(CRC_Fib$Celltype))
CRC_Fib$orig.ident <- as.factor(as.character(CRC_Fib$orig.ident))
saveRDS(CRC_Fib, file = "CRC_HRA000979_Fib.rds")

#CAFs
CRC_Fib_T <- subset(CRC_Fib, subset = Sample %in% c("Tumor"))
CRC_Fib_T$Sample <- as.factor(as.character(CRC_Fib_T$Sample))
CRC_Fib_T$orig.ident <- as.factor(as.character(CRC_Fib_T$orig.ident))
saveRDS(CRC_Fib_T, file = "CRC_HRA000979_Fib_T.rds")
Fib <- CRC_Fib_T
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
saveRDS(Fib, file = "CRC_HRA000979_Fib_T_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4E.
library(clustree)
library(patchwork)
pdf(file = "clustree_CRC_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
Fib <- FindClusters(Fib, resolution = 0.2)
saveRDS(Fib, file = "CRC_HRA000979_Fib_T_res0.2.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="CRC_Fib_T_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="CRC_Fib_T_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="CRC_Fib_T_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="CRC_Fib_T_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5", "MYH11")
pdf(file="CRC_Fib_T_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(2)
myCAFs = c(0, 3)
MYH11_CAFs = c(4)
iCAFs = c(1)
current.cluster.ids <- c(apCAFs,
                         myCAFs,
                         MYH11_CAFs,
                         iCAFs)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
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
saveRDS(Fib, file = "CRC_HRA000979_Fib_T_res0.2_annotation.rds")

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "MYH11+ CAFs" = "MYH11+ CAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "CRC_HRA000979_Fib_T_res0.2_annotation.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5F: left.
pdf(file="CRC_HRA000979_Fib_T_umap_celltype.pdf",width=5.7, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5F: right.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "MYH11+ CAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "CXCL12", "CXCL1", "CCL2", "MYH11", "RGS5", "ACTA2", "HLA-DRB1", "HLA-DRA", "CD74")
pdf(file="CRC_HRA000979_Fib_T_F_features_res0.2_annotation.pdf",width=8, height=5.5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
