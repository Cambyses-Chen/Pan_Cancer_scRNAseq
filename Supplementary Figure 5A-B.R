###codes for Supplementary Figure 5A-B.

#1 Supplementary Figure 5A.
#NPC.
###Read and preprocess the original count matrix.
library(Seurat)
library(ggplot2)
library(scales)
library(ggsci)
##1
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N52_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder),
                     project = "NPC" )
})
sceList
N52 <- sceList[[1]]
N52
table(N52$orig.ident)
head(colnames(N52))
##2
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N53_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N53 <- sceList[[1]]
N53
table(N53$orig.ident)
head(colnames(N53))

##3
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N54_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N54 <- sceList[[1]]
N54
table(N54$orig.ident)
head(colnames(N54))

##4
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N55_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N55 <- sceList[[1]]
N55
table(N55$orig.ident)
head(colnames(N55))

##5
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N57_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N57 <- sceList[[1]]
N57
table(N57$orig.ident)
head(colnames(N57))

##6
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N61_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N61 <- sceList[[1]]
N61
table(N61$orig.ident)
head(colnames(N61))

##7
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/N68_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
N68 <- sceList[[1]]
N68
table(N68$orig.ident)
head(colnames(N68))

##8
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC33_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
###Create one Seurat object
NPC33 <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
               project = "NPC")
NPC33
table(NPC33$orig.ident)
head(colnames(NPC33))
tail(colnames(NPC33))

##9
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC36_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
###Create one Seurat object
NPC36 <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
               project = "NPC")
NPC36
table(NPC36$orig.ident)
head(colnames(NPC36))
tail(colnames(NPC36))

##10
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC46_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
###Create one Seurat object
NPC46 <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
               project = "NPC")
NPC46
table(NPC46$orig.ident)
head(colnames(NPC46))
tail(colnames(NPC46))

##11
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC47_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
###Create one Seurat object
NPC47 <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
               project = "NPC")
NPC47
table(NPC47$orig.ident)
head(colnames(NPC47))
tail(colnames(NPC47))

##12
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC49_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
###Create one Seurat object
NPC49 <- merge(sceList[[1]], 
               y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
               project = "NPC")
NPC49
table(NPC49$orig.ident)
head(colnames(NPC49))
tail(colnames(NPC49))

##13
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC50_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC50 <- sceList[[1]]
NPC50
table(NPC50$orig.ident)
head(colnames(NPC50))

##14
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC51_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC51 <- sceList[[1]]
NPC51
table(NPC51$orig.ident)
head(colnames(NPC51))

##15
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC58_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC58 <- sceList[[1]]
NPC58
table(NPC58$orig.ident)
head(colnames(NPC58))

##16
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC59_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC59 <- sceList[[1]]
NPC59
table(NPC59$orig.ident)
head(colnames(NPC59))

##17
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC60_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC60 <- sceList[[1]]
NPC60
table(NPC60$orig.ident)
head(colnames(NPC60))

##18
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC62_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC62 <- sceList[[1]]
NPC62
table(NPC62$orig.ident)
head(colnames(NPC62))

##19
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC63_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC63 <- sceList[[1]]
NPC63
table(NPC63$orig.ident)
head(colnames(NPC63))

##20
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC64_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC64 <- sceList[[1]]
NPC64
table(NPC64$orig.ident)
head(colnames(NPC64))

##21
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC65_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC65 <- sceList[[1]]
NPC65
table(NPC65$orig.ident)
head(colnames(NPC65))

##22
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC66_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC66 <- sceList[[1]]
NPC66
table(NPC66$orig.ident)
head(colnames(NPC66))

##23
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087/NPC69_10X")
folders=list.files('./','^HRR')
folders
sceList = lapply(folders,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = "NPC" )
})
sceList
NPC69 <- sceList[[1]]
NPC69
table(NPC69$orig.ident)
head(colnames(NPC69))

###annotation
N52$Sample = "Normal"
N53$Sample = "Normal"
N54$Sample = "Normal"
N55$Sample = "Normal"
N57$Sample = "Normal"
N61$Sample = "Normal"
N68$Sample = "Normal"
NPC33$Sample = "Tumor"
NPC36$Sample = "Tumor"
NPC46$Sample = "Tumor"
NPC47$Sample = "Tumor"
NPC49$Sample = "Tumor"
NPC50$Sample = "Tumor"
NPC51$Sample = "Tumor"
NPC58$Sample = "Tumor"
NPC59$Sample = "Tumor"
NPC60$Sample = "Tumor"
NPC62$Sample = "Tumor"
NPC63$Sample = "Tumor"
NPC64$Sample = "Tumor"
NPC65$Sample = "Tumor"
NPC66$Sample = "Tumor"
NPC69$Sample = "Tumor"

###Create one Seurat object
NPC  <- merge(N52, y = c(N53, N54, N55, N57, N61, N68, NPC33, NPC36, NPC46, NPC47, NPC49, NPC50, NPC51, NPC58, NPC59, NPC60, NPC62, NPC63, NPC64, NPC65, NPC66, NPC69), project = "NPC")
NPC
head(colnames(NPC))
tail(colnames(NPC))
table(NPC$orig.ident)
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
saveRDS(NPC, file = "NPC_seurat.rds")
dim(NPC)

###Preliminary filtering.
selected_c <- WhichCells(NPC, expression = nFeature_RNA > 200)
selected_f <- rownames(NPC)[Matrix::rowSums(NPC) > 3]
NPC <- subset(NPC, features = selected_f, cells = selected_c)
dim(NPC)

###QC metrics and filtering.
NPC[["percent.mt"]] <- PercentageFeatureSet(NPC, pattern = "^MT-")
pdf(file="NPC.featureViolin.pdf",width=18,height=6)
VlnPlot(object = NPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pdf(file="NPC.FeatureScatter.pdf",width=12,height=6)
plot1 <- FeatureScatter(NPC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(NPC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
NPC <- subset(NPC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
dim(NPC)

###Visualization.
#VlnPlot for Supplementary Figure 1E.
Idents(NPC) <- NPC$orig.ident
pdf(file="NPC.featureViolin_QC.pdf",width=18,height=6)           
VlnPlot(object = NPC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

###Normalizing the data.
NPC <- NormalizeData(NPC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC <- FindVariableFeatures(NPC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC)
NPC <- ScaleData(NPC, features = all.genes)
NPC <- ScaleData(NPC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC <- RunPCA(NPC, features = VariableFeatures(object = NPC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC <- RunHarmony(NPC, group.by.vars = "orig.ident")
NPC <- RunUMAP(NPC, reduction = "harmony", dims = 1:20)
NPC <- RunTSNE(NPC, reduction = "harmony", dims = 1:20)
names(NPC@reductions)

###Cluster the cells.
NPC <- FindNeighbors(NPC, reduction = "harmony", dims = 1:20)
saveRDS(NPC, file = "NPC_before_resolution.rds")

#Determine resolution.
test <- NPC
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 3C.
library(clustree)
library(patchwork)
pdf(file = "clustree_NPC.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.2.
NPC <- FindClusters(NPC, resolution = 0.2)
saveRDS(NPC, file = "NPC_resolution_0.2.rds")

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
pdf(file="NPC_marker.pdf",width=44, height=33)
VlnPlot(NPC, features = celltype_marker, pt.size = 0, ncol = 2)
dev.off()

##2
NK_marker=c("FGFBP2", "KLRD1", "CX3CR1", "GNLY", "NKG7", "CD3D", "GZMA", "GZMB")
pdf(file="NPC_NK.pdf",width=44, height=14)
VlnPlot(NPC, features = NK_marker, pt.size = 0, ncol = 2)
dev.off()

##3
Epithelial=c("EPCAM", "KRT19", "PROM1", "ALDH1A1", "CD24")
pdf(file="NPC_Epithelial.pdf",width=44, height=14)
VlnPlot(NPC, features = Epithelial, pt.size = 0, ncol = 2)
dev.off()

##4
Plasma_cells=c("IGHG1", "MZB1", "SDC1", "CD79A", "JCHAIN")
pdf(file="NPC_Plasma_cells.pdf",width=44, height=14)
VlnPlot(NPC, features = Plasma_cells, pt.size = 0, ncol = 2)
dev.off()

##5
B_marker=c("CD19", "CD79A", "CD79A", "MS4A1", "CD27", "CD1D", "CD38")
pdf(file="NPC_B_marker.pdf",width=44, height=19)
VlnPlot(NPC, features = B_marker, pt.size = 0, ncol = 2)
dev.off()

##6
Fibroblasts=c("FGF7", "MME", "COL3A1", "PDGFRA", "RGS5", "COL1A1")
pdf(file="NPC_Fibroblasts.pdf",width=44, height=14)
VlnPlot(NPC, features = Fibroblasts, pt.size = 0, ncol = 2)
dev.off()

##7
Endothelial=c("PECAM1", "VWF")
pdf(file="NPC_Endothelial.pdf",width=44, height=5)
VlnPlot(NPC, features = Endothelial, pt.size = 0, ncol = 2)
dev.off()

##8
T_marker=c("CD3D", "CD3E", "CD8A", "CD4", "SELL", "IL2RA")
pdf(file="NPC_T_marker.pdf",width=44, height=14)
VlnPlot(NPC, features = T_marker, pt.size = 0, ncol = 2)
dev.off()

##9
Mast=c("KIT", "TPSAB1")
pdf(file="NPC_Mast.pdf",width=44, height=5)
VlnPlot(NPC, features = Mast, pt.size = 0, ncol = 2)
dev.off()

##10
myeloid_cells=c("LYZ", "CD163", "AIF1")
pdf(file="NPC_myeloid.pdf",width=44, height=9.5)
VlnPlot(NPC, features = myeloid_cells, pt.size = 0, ncol = 2)
dev.off()

##11
Monocyte_macrophage=c("CD68", "CD163", "CD14", "MRC1")
pdf(file="NPC_Monocyte_macrophage.pdf",width=44, height=9.5)
VlnPlot(NPC, features = Monocyte_macrophage, pt.size = 0, ncol = 2)
dev.off()

##12
pDC = c("PTGDS", "GZMB", "IGJ", "TCL1A", "LILRA4", "PLAC8", "ITM2C", "IRF7")
pdf(file="NPC_pDC.pdf",width=44, height=20)
VlnPlot(NPC, features = pDC, pt.size = 0, ncol = 2)
dev.off()

##13
neutrophil = c("CEACAM8", "FCGR3B")
pdf(file="NPC_neutrophil.pdf",width=44, height=5)
VlnPlot(NPC, features = neutrophil, pt.size = 0, ncol = 2)
dev.off()

##14
melanocyte=c("SOX10", "MITF", "DCT", "MLANA", "PMEL", "TYR", "TYRP1")
pdf(file="NPC_melanocyte.pdf",width=44, height=30)
VlnPlot(NPC, features = melanocyte, pt.size = 0, ncol = 2)
dev.off()

##15
clear_cell=c("CA9", "GPX3", "GATM")
pdf(file="NPC_clear_cell.pdf",width=44, height=10)
VlnPlot(NPC, features = clear_cell, pt.size = 0, ncol = 2)
dev.off()

###Assigning cell type identity to clusters.
Fibroblasts = c(8)
Myeloid_cells = c(5)
Mast_cells = c(9)
pDCs = c(11)
B_cells = c(0, 13)
Plasma_cells = c(7)
T_cells = c(1, 2, 4, 6)
Endothelial_cells = c(10)
Epithelial_cells = c(3)
unknown = c(12)
current.cluster.ids <- c(Fibroblasts,
                         Myeloid_cells,
                         Mast_cells,
                         pDCs,
                         B_cells,
                         Plasma_cells,
                         T_cells,
                         Endothelial_cells,
                         Epithelial_cells,
                         unknown)
new.cluster.ids <- c(rep("Fibroblasts",length(Fibroblasts)),
                     rep("Myeloid_cells",length(Myeloid_cells)),
                     rep("Mast_cells",length(Mast_cells)),
                     rep("pDCs",length(pDCs)),
                     rep("B_cells",length(B_cells)),
                     rep("Plasma_cells",length(Plasma_cells)),
                     rep("T_cells",length(T_cells)),
                     rep("Endothelial_cells",length(Endothelial_cells)),
                     rep("Epithelial_cells",length(Epithelial_cells)),
                     rep("unknown",length(unknown))
)

NPC@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(NPC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(NPC$Celltype)
table(NPC$seurat_clusters)

###Rename.
Idents(NPC) <- NPC$Celltype
NPC <- RenameIdents(NPC, "T_cells" = "T cells", "Plasma_cells" = "Plasma cells", "B_cells" = "B cells", "Epithelial_cells" = "Epithelial cells", "Myeloid_cells" = "Myeloid cells", "Mast_cells" = "Mast cells", "Endothelial_cells" = "Endothelial cells")
table(Idents(NPC))
NPC$Celltype <- Idents(NPC)
table(NPC$Celltype)
saveRDS(NPC, file = "NPC_resolution0.2_annotation.rds")

###Exclude unknown cell subpopulations.
NPC <- subset(NPC, subset = Celltype %in% c("Fibroblasts", "T cells", "Plasma cells", "B cells", "Epithelial cells", "Myeloid cells", "Mast cells", "Endothelial cells", "pDCs"))
table(NPC$Celltype)
NPC$Celltype <- as.factor(as.character(NPC$Celltype))
NPC$orig.ident <- as.factor(as.character(NPC$orig.ident))
NPC$Sample <- as.factor(as.character(NPC$Sample))

###Ranking.
Idents(NPC) <- NPC$Celltype
NPC <- RenameIdents(NPC, "T cells" = "T cells", "B cells" = "B cells", "Plasma cells" = "Plasma cells", "pDCs" = "pDCs", "Mast cells" = "Mast cells", "Myeloid cells" = "Myeloid cells", "Fibroblasts" = "Fibroblasts", "Endothelial cells" = "Endothelial cells", "Epithelial cells" = "Epithelial cells")
table(Idents(NPC))
NPC$Celltype <- Idents(NPC)
table(NPC$Celltype)
saveRDS(NPC, file = "NPC_resolution0.2_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5A: left.
pdf(file="NPC_umap_Celltype.pdf",width=5.8, height=4)
DimPlot(NPC, reduction = "umap",label=T, group.by = "Celltype", repel = T, raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5A: right.
Idents(NPC) <- NPC$Celltype
My_levels <- c( "T cells", "B cells", "Plasma cells", "pDCs", "Mast cells", "Myeloid cells", "Fibroblasts", "Endothelial cells", "Epithelial cells")
Idents(NPC) <- factor(Idents(NPC), levels= My_levels)
features = c("KRT19", "EPCAM", "VWF", "PECAM1", "COL3A1", "COL1A1" , "AIF1", "LYZ", "TPSAB1", "KIT", "TCL1A", "LILRA4", "IGHG1", "JCHAIN", "CD79A", "CD19", "MS4A1", "CD3E", "CD3D")
pdf(file="NPC_features.pdf",width=8, height=9)
DotPlot(NPC, features = features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

#DimPlot for Supplementary Figure 2C.
pdf(file="NPC_umap_orig.ident.pdf",width=6, height=4)
DimPlot(NPC, reduction = "umap",label=T, group.by = "orig.ident", raster=FALSE) + ggtitle("") + scale_color_igv()
dev.off()

#2 Supplementary Figure 5B.
###Cluster the Fibroblasts in tumor.
NPC_Fib <- subset(NPC, subset = Celltype %in% c("Fibroblasts"))
NPC_Fib$seurat_clusters <- as.factor(as.character(NPC_Fib$seurat_clusters))
NPC_Fib$Celltype <- as.factor(as.character(NPC_Fib$Celltype))
NPC_Fib$orig.ident <- as.factor(as.character(NPC_Fib$orig.ident))
NPC_Fib$Sample <- as.factor(as.character(NPC_Fib$Sample))
saveRDS(NPC_Fib, file = "NPC_Fib.rds")

#CAFs.
NPC_Fib_T <- subset(NPC_Fib, subset = Sample %in% c("Tumor"))
NPC_Fib_T$seurat_clusters <- as.factor(as.character(NPC_Fib_T$seurat_clusters))
NPC_Fib_T$Celltype <- as.factor(as.character(NPC_Fib_T$Celltype))
NPC_Fib_T$orig.ident <- as.factor(as.character(NPC_Fib_T$orig.ident))
NPC_Fib_T$Sample <- as.factor(as.character(NPC_Fib_T$Sample))
saveRDS(NPC_Fib_T, file = "NPC_Fib_T.rds")
Fib <- NPC_Fib_T
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
saveRDS(Fib, file = "NPC_Fib_T_before.rds")

#Determine resolution.
test <- Fib
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  test <- FindClusters(test, resolution = res)
}

###Visualization.
#clustering tree for Supplementary Figure 4C.
library(clustree)
library(patchwork)
pdf(file = "clustree_NPC_Fib.pdf", width = 8, height = 10)
clustree(test, prefix = 'RNA_snn_res.') + coord_flip()
dev.off()

#resolution = 0.1.
Fib <- FindClusters(Fib, resolution = 0.1)
saveRDS(Fib, file = "NPC_Fib_T_res0.1.rds")

###classic celltype marker expression across different clusters.
apCAFs = c("HLA-DRA", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74")
pdf(file="NPC_Fib_T_apCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = apCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

vCAFs = c("VEGFA", "ACTA2", "LRRC15", "NID2")
pdf(file="NPC_Fib_T_vCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = vCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

iCAFs = c("CD34", "CXCL1", "CXCL12", "CXCL13",  "IL6", "CXCL14", "CCL2", "PDGFRA", "HAS1")
pdf(file="NPC_Fib_T_iCAFs.pdf",width=8, height=5)
DotPlot(Fib, features = iCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

eCAFs = c("MMP14", "LOXL2", "POSTN")
pdf(file="NPC_Fib_T_eCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = eCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

myCAFs = c("ACTA2", "CCN2", "POSTN", "RGS5")
pdf(file="NPC_Fib_T_myCAFs.pdf",width=8, height=4)
DotPlot(Fib, features = myCAFs, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()

###Assigning cell type identity to clusters.
apCAFs = c(5)
iCAFs = c(4)
myCAFs = c(1)
eCAFs = c(0, 3, 6)
unknown_CAFs = c(2)
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
saveRDS(Fib, file = "NPC_Fib_T_res0.1_annotation.rds")

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
saveRDS(Fib, file = "NPC_Fib_T_res0.1_annotation_sub.rds")

###Visualization.
#UMAP plot for Supplementary Figure 5B: left.
pdf(file="NPC_T_Fib_umap_celltype.pdf",width=5.3, height=4)
DimPlot(Fib, reduction = "umap",label=T, group.by = "celltype", repel = T) + ggtitle("") + scale_color_igv()
dev.off()

#Bubble heatmap for Supplementary Figure 5B: right.
Idents(Fib) <- 'celltype'
My_levels <- c( "apCAFs", "myCAFs", "eCAFs", "iCAFs")
Idents(Fib) <- factor(Idents(Fib), levels= My_levels)
F_features = c("CXCL14", "CXCL12", "CD34", "POSTN", "MMP14", "RGS5", "ACTA2", "HLA-DQB1", "HLA-DPA1", "HLA-DRA", "CD74")
pdf(file="NPC_Fib_T_F_features_annotation.pdf",width=8, height=6)
DotPlot(Fib, features = F_features, dot.scale = 12)+RotatedAxis()+coord_flip()+xlab('Gene')+ylab('Cluster')+FontSize(x.text = 20, y.text = 20, x.title = 20, y.title = 20, main = 20) + scale_colour_gradient2(low = muted("blue"), mid = "white", high = muted("red"))
dev.off()
