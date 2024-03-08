###codes for Supplementary Figure 5G-H.

#1 Supplementary Figure 5G.
#BRCA.
###Read and preprocess the original count matrix.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
BC <- Read10X(data.dir = "./BrCa_Atlas_Count_out/")
BC <- CreateSeuratObject(counts = BC, project = "BC")
BC
head(BC@meta.data)
table(BC$orig.ident)

###Split the Seurat object.
scRNAlist <- SplitObject(BC, split.by = "orig.ident")
scRNAlist

###annotation.
scRNAlist$CID3586$Sample <- "HER2"
scRNAlist$CID3921$Sample <- "HER2"
scRNAlist$CID45171$Sample <- "HER2"
scRNAlist$CID3838$Sample <- "HER2"
scRNAlist$CID4066$Sample <- "HER2"
scRNAlist$CID44041$Sample <- "TNBC"
scRNAlist$CID4465$Sample <- "TNBC"
scRNAlist$CID4495$Sample <- "TNBC"
scRNAlist$CID44971$Sample <- "TNBC"
scRNAlist$CID44991$Sample <- "TNBC"
scRNAlist$CID4513$Sample <- "TNBC"
scRNAlist$CID4515$Sample <- "TNBC"
scRNAlist$CID4523$Sample <- "TNBC"
scRNAlist$CID3946$Sample <- "TNBC"
scRNAlist$CID3963$Sample <- "TNBC"
scRNAlist$CID4461$Sample <- "ER"
scRNAlist$CID4463$Sample <- "ER"
scRNAlist$CID4471$Sample <- "ER"
scRNAlist$CID4530N$Sample <- "ER"
scRNAlist$CID4535$Sample <- "ER"
scRNAlist$CID4040$Sample <- "ER"
scRNAlist$CID3941$Sample <- "ER"
scRNAlist$CID3948$Sample <- "ER"
scRNAlist$CID4067$Sample <- "ER"
scRNAlist$CID4290A$Sample <- "ER"
scRNAlist$CID4398$Sample <- "ER"

###Create one Seurat object
BC <- merge(x = scRNAlist$CID3586, y = c(scRNAlist$CID3921,
                                         scRNAlist$CID45171,
                                         scRNAlist$CID3838,
                                         scRNAlist$CID4066,
                                         scRNAlist$CID44041,
                                         scRNAlist$CID4465,
                                         scRNAlist$CID4495,
                                         scRNAlist$CID44971,
                                         scRNAlist$CID44991,
                                         scRNAlist$CID4513,
                                         scRNAlist$CID4515,
                                         scRNAlist$CID4523,
                                         scRNAlist$CID3946,
                                         scRNAlist$CID3963,
                                         scRNAlist$CID4461,
                                         scRNAlist$CID4463,
                                         scRNAlist$CID4471,
                                         scRNAlist$CID4530N,
                                         scRNAlist$CID4535,
                                         scRNAlist$CID4040,
                                         scRNAlist$CID3941,
                                         scRNAlist$CID3948,
                                         scRNAlist$CID4067,
                                         scRNAlist$CID4290A,
                                         scRNAlist$CID4398), project = "BC")
BC
head(colnames(BC))
tail(colnames(BC))
table(BC$orig.ident)
table(BC$Sample)
saveRDS(BC, file = "BC_GSE176078_Merge.rds")

###Preliminary filtering.
selected_c <- WhichCells(BC, expression = nFeature_RNA > 200)
selected_f <- rownames(BC)[Matrix::rowSums(BC) > 3]
BC <- subset(BC, features = selected_f, cells = selected_c)
dim(BC)

###QC metrics and filtering.
BC[["percent.mt"]] <- PercentageFeatureSet(BC, pattern = "^MT-")
Idents(BC) <- BC$orig.ident
pdf(file="BC.featureViolin.pdf",width=18,height=6)           
VlnPlot(object = BC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="BC.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(BC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
BC <- subset(BC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(BC)

###Visualization.
#VlnPlot for Supplementary Figure 1G.
Idents(BC) <- BC$orig.ident
pdf(file="BC.featureViolin_QC.pdf",width=16,height=6)           
VlnPlot(object = BC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
BC <- NormalizeData(BC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
BC <- FindVariableFeatures(BC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(BC)
BC <- ScaleData(BC, features = all.genes)
BC <- ScaleData(BC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
BC <- RunPCA(BC, features = VariableFeatures(object = BC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
BC <- RunHarmony(BC, group.by.vars = "orig.ident")
BC <- RunUMAP(BC, reduction = "harmony", dims = 1:20)
BC <- RunTSNE(BC, reduction = "harmony", dims = 1:20)
names(BC@reductions)

###Cluster the cells.
BC <- FindNeighbors(BC, reduction = "harmony", dims = 1:20)
saveRDS(BC, file = "BC_GSE176078_before_resolution.rds")

#Determine resolution.
test <- BC
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3F.
library(clustree)
library(patchwork)
pdf(file = "clustree_BC.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
BC <- FindClusters(BC, resolution = 0.1)
saveRDS(BC, file = "BC_GSE176078_resolution_0.1.rds")

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
pdf(file="BC_marker.pdf",width=44, height=33)
VlnPlot(BC, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="BC_NK.pdf",width=44, height=14)
VlnPlot(BC, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="BC_Epithelial.pdf",width=44, height=14)
VlnPlot(BC, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="BC_Plasma_cells.pdf",width=44, height=14)
VlnPlot(BC, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="BC_B_marker.pdf",width=44, height=19)
VlnPlot(BC, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="BC_Fibroblasts.pdf",width=44, height=14)
VlnPlot(BC, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="BC_Endothelial.pdf",width=44, height=5)
VlnPlot(BC, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="BC_T_marker.pdf",width=44, height=14)
VlnPlot(BC, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="BC_Mast.pdf",width=44, height=5)
VlnPlot(BC, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="BC_myeloid.pdf",width=44, height=9.5)
VlnPlot(BC, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="BC_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(BC, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="BC_pDC.pdf",width=44, height=20)
VlnPlot(BC, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="BC_neutrophil.pdf",width=44, height=5)
VlnPlot(BC, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="BC_melanocyte.pdf",width=44, height=30)
VlnPlot(BC, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="BC_clear_cell.pdf",width=44, height=10)
VlnPlot(BC, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(5, 6)
Myeloid_cells = c(3)
B_cells = c(8)
Plasma_cells = c(7)
T_cells = c(0, 1, 9)
Epithelial_cells = c(2, 11)
Endothelial_cells = c(4)
unknown = c(10)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Epithelial_cells,
                         Endothelial_cells,
                         unknown)

new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("unknown",length(unknown))
)

BC@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(BC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(BC@meta.data$Celltype)
table(BC$seurat_clusters)

###Rename.
Idents(BC) <- BC$Celltype
BC <- RenameIdents(BC, "T_cells" = "T cells", "Endothelial_cells" = "Endothelial cells", "Plasma_cells" = "Plasma cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells")
table(Idents(BC))
BC$Celltype <- Idents(BC)
table(BC$Celltype)
saveRDS(BC, file = "BC_GSE176078_resolution0.1_annotation.rds")

###Exclude unknown cell subpopulations.
BC <- subset(BC, subset = Celltype %in% c("Fibroblasts", "Myeloid cells", "B cells", "Plasma cells", "T cells", "Epithelial cells", "Endothelial cells"))
BC$seurat_clusters <- as.factor(as.character(BC$seurat_clusters))
BC$Celltype <- as.factor(as.character(BC$Celltype))
BC$orig.ident <- as.factor(as.character(BC$orig.ident))
BC$Sample <- as.factor(as.character(BC$Sample))

###Ranking.
Idents(BC) <- BC$Celltype
BC <- RenameIdents(BC, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Epithelial cells" = "Epithelial cells")
table(Idents(BC))
BC$Celltype <- Idents(BC)
table(BC$Celltype)
saveRDS(BC, file = "BC_GSE176078_resolution0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5G: left.
pdf(file="BC_GSE176078_umap_Celltype.pdf",width=5.7, height=4)
DimPlot(BC, reduction = "umap",label=T, group.by = "Celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5G: right.
Idents(BC) <- 'Celltype'
My_levels <- c("T cells", "B cells", "Plasma cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Epithelial cells")
Idents(BC) <- factor(Idents(BC), levels= My_levels)
features = c("KRT19", "EPCAM", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="BC_GSE176078_features.pdf",width=8, height=8)
DotPlot(BC, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2F.
pdf(file="BC_umap_orig.ident.pdf",width=6.3, height=4)
DimPlot(BC, reduction = "umap",label=T, group.by = "orig.ident") + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5H.
###Cluster the Fibroblasts
BC_Fib <- subset(BC, subset = Celltype %in% c("Fibroblasts"))
table(BC_Fib$Celltype)
BC_Fib$seurat_clusters <- as.factor(as.character(BC_Fib$seurat_clusters))
BC_Fib$Celltype <- as.factor(as.character(BC_Fib$Celltype))
saveRDS(BC_Fib, file = "BC_GSE176078_Fib.rds")
Fib <- BC_Fib
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
saveRDS(Fib, file = "BC_GSE176078_Fib_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4F.
library(clustree)
library(patchwork)
pdf(file = "clustree_BC_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
Fib <- FindClusters(Fib, resolution = 0.2)
saveRDS(Fib, file = "BC_GSE176078_Fib_res0.2.rds")

###classic celltype marker in clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="BC_GSE176078_Fib_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="BC_GSE176078_Fib_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="BC_GSE176078_Fib_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="BC_GSE176078_Fib_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5", "MYH11")
pdf(file="BC_GSE176078_Fib_myCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(5)
eCAFs = c(1)
myCAFs = c(0, 2)
iCAFs = c(3, 4)
current.cluster.ids <- c(apCAFs,
                         eCAFs,
                         myCAFs,
                         iCAFs)
new.cluster.ids <- c(rep("apCAFs",length(apCAFs)),
                     rep("eCAFs",length(eCAFs)),
                     rep("myCAFs",length(myCAFs)),
                     rep("iCAFs",length(iCAFs))
)
Fib@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(Fib@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(Fib$celltype)
table(Fib$seurat_clusters)
saveRDS(Fib, file = "BC_GSE176078_Fib_res0.2_annotation.rds")

###Ranking.
Idents(Fib) <- Fib$celltype
Fib <- RenameIdents(Fib, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "iCAFs" = "iCAFs")
table(Idents(Fib))
Fib$celltype <- Idents(Fib)
table(Fib$celltype)
saveRDS(Fib, file = "BC_GSE176078_Fib_res0.2_annotation2.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5H: left.
pdf(file="BC_GSE176078_Fib_umap_celltype.pdf",width=5.3, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5H: right.
Idents(Fib) <- Fib$celltype
F_features = c("CXCL14", "CXCL12", "CXCL1", "CD34", "POSTN", "MMP14", "RGS5", "ACTA2", "HLA-DRB1", "HLA-DRA", "CD74")
pdf(file="BC_GSE176078_Fib_F_features_res0.2_annotation.pdf",width=8, height=5.5)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
