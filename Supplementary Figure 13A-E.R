###codes for Supplementary Figure 13A-E.

library(semla) 
library(tibble) 
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#1 Supplementary Figure 13A, B, C.
#GSM7089855 slice.
###Load CRC Spatial transcriptome (ST) data.
#Read the file path.
setwd("E:/Raw Data/ST_data/CRC")
samples <- Sys.glob("./GSM7089855/Spatial/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("./GSM7089855/Spatial/tissue_hires_image.png")
spotfiles <- Sys.glob("./GSM7089855/Spatial/tissue_positions_list.csv")
json <- Sys.glob("./GSM7089855/Spatial/scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, 
                    section_id = c("GSM7089855"))

###Create ST Seurat object.
GSM7089855 <- ReadVisiumData(infoTable)
saveRDS(GSM7089855, file = "GSM7089855_CRC_semla.rds")

###Load the annotated scRNA-seq data for reference.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
SC_CRC <- readRDS(file = "CRC_Tumor_merge.rds") #CRC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of CRC tumors. The specific code is in the file named CRC_Tumor_merge.R.
table(SC_CRC$celltype)

###Normalize data and find variable features for ST data.
GSM7089855 <- GSM7089855 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000)

###Normalize scRNA-seq data and run vanilla analysis to create UMAP embedding.
SC_CRC <- SC_CRC |>
  NormalizeData() |>
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA() |> 
  RunUMAP(reduction = "pca", dims = 1:30)

###Rerun FindVariableFeatures to increase the number before cell type deconvolution.
SC_CRC <- SC_CRC |> 
  FindVariableFeatures(nfeatures = 10000)

###Visualize the available cell types in our UMAP embedding of the cells.
DimPlot(SC_CRC, group.by = "celltype")

###Run NNLS.(a quick method based on Non-Negative Least Squares to infer cell type proportions directly from Visium spot expression profiles.)
DefaultAssay(GSM7089855) <- "Spatial"
ti <- Sys.time()
GSM7089855 <- RunNNLS(object = GSM7089855, 
                      singlecell_object = SC_CRC, 
                      groups = "celltype")
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))

###Check available cell types.
rownames(GSM7089855)

###Selected cell types.
DefaultAssay(GSM7089855) <- "celltypeprops"
selected_celltypes <- rownames(GSM7089855)

###Loading H&E images.
setwd("E:/Raw Data/ST_data/CRC")
GSM7089855 <- LoadImages(GSM7089855, image_height = 1e3)
saveRDS(GSM7089855, file = "GSM7089855_CRC_semla2.rds")

###Visualization.
#Visualize single cell type.
setwd("E:/Raw Data/ST_data/CRC")
#MapFeatures plot for Supplementary Figure 13A: left.
pdf(file="CRC_ST_MapFeatures_GSM7089855_Tumor_cells.pdf",width=4,height=5)
MapFeatures(GSM7089855, pt_size = 1.7,
            features = selected_celltypes[13], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#MapFeatures plot for Supplementary Figure 13A: right.
pdf(file="CRC_ST_MapFeatures_GSM7089855_apCAFs.pdf",width=4,height=5)
MapFeatures(GSM7089855, pt_size = 1.7,
            features = selected_celltypes[15], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#Visualize multiple cell types.
#Loading H&E images.
setwd("E:/Raw Data/ST_data/CRC")
GSM7089855 <- GSM7089855 |> 
  LoadImages()

#MapMultipleFeatures plot for Supplementary Figure 13B.
pdf(file="CRC_ST_MapMultipleFeatures_GSM7089855_Tumor_cells_apCAFs.pdf",width=6,height=4)
MapMultipleFeatures(GSM7089855,
                    image_use = "raw",
                    section_number = 1,
                    label_by = "section_id",
                    ncol = 1,
                    pt_size = 1.7, max_cutoff = 0.99,
                    override_plot_dims = TRUE, 
                    colors = c("#88CCEE", "#CC6677"),
                    features = selected_celltypes[c(13, 15)])
dev.off()

###Obtain the selected cell type proportions for each spot in the tissue slice.
t <- FetchData(GSM7089855, selected_celltypes[c(13, 15)])
names(t) <- c("Tumor_cells", "apCAFs")#rename cell type name.

###Correlation analysis.
data <- t
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor_cells)

###Visualization.
#Scatter plot for Supplementary Figure 13C.
setwd("E:/Raw Data/ST_data/CRC")
pdf(file = "CRC_ST_GSM7089855_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor_cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 13D, E.
###Load CRC Spatial transcriptome (ST) data.
#GSM7089855 slice.
setwd("E:/Raw Data/ST_data/CRC/GSM7089855")
t <- Read10X_h5("filtered_feature_bc_matrix.h5")
GSM7089855 <- Seurat::CreateSeuratObject(counts = t, project = 'GSM7089855', assay = 'Spatial')
GSM7089855
GSM7089855$slice <- "GSM7089855"
S7089855 <- Seurat::Read10X_Image(image.dir = './Spatial/')
Seurat::DefaultAssay(object = S7089855) <- 'Spatial'
S7089855 <- S7089855[colnames(x = GSM7089855)]
GSM7089855[['GSM7089855']] <- S7089855
GSM7089855
setwd("E:/Raw Data/ST_data/CRC")
saveRDS(GSM7089855, file = "CRC_ST_Seurat_GSM7089855.rds")

###Apply sctransform normalization.
GSM7089855  <- SCTransform(GSM7089855 , assay = "Spatial", verbose = FALSE)
saveRDS(GSM7089855, file = "CRC_ST_Seurat_GSM7089855_sctransform.rds")

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/CRC")
gs = lapply(readLines("CRC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
GSM7089855 <- AddModuleScore(
  object = GSM7089855,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
head(GSM7089855@meta.data)
colnames(GSM7089855@meta.data)[7] <- 'apCAFs signature'
colnames(GSM7089855@meta.data)[8] <- 'CD4+ effector T cells signature'

##Visualization.
setwd("E:/Raw Data/ST_data/CRC")
#SpatialFeaturePlot for Supplementary Figure 13D: left.
pdf(file="CRC_ST_SpatialFeaturePlot_GSM7089855_CD4_effector_T.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(GSM7089855, features = "CD4+ effector T cells signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "GSM7089855")
dev.off()

#SpatialFeaturePlot for Supplementary Figure 13D: right.
pdf(file="CRC_ST_SpatialFeaturePlot_GSM7089855_apCAFs.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(GSM7089855, features = "apCAFs signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "GSM7089855")
dev.off()

###Obtain module scores for cell type signatures for each spot in the tissue slice.
m <- GSM7089855@meta.data
df_GSM7089855 <- m[, c(4, 7, 8)]

###Rename.
colnames(df_GSM7089855) <- c("slice", "apCAFs", "CD4T")

###Correlation analysis.
#Shapiro-Wilk test for normality.
shapiro.test(df_GSM7089855$CD4T)
shapiro.test(df_GSM7089855$apCAFs)

###Visualization.
setwd("E:/Raw Data/ST_data/CRC")
#scatter plot for Supplementary Figure 13E.
pdf(file = "CRC_ST_GSM7089855_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_GSM7089855, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
