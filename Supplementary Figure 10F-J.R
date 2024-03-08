###codes for Supplementary Figure 10F-J.

library(semla) 
library(tibble) 
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#1 Supplementary Figure 10F, G, H.
#P210325T3 slice.
###Load HNSCC Spatial transcriptome (ST) data.
#Read the file path.
setwd("E:/Raw Data/ST_data/HNSCC")
samples <- Sys.glob("./P210325T3/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("./P210325T3/tissue_hires_image.png")
spotfiles <- Sys.glob("./P210325T3/tissue_positions_list.csv")
json <- Sys.glob("./P210325T3/scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, 
                    section_id = c("P210325T3"))

###Create ST Seurat object.
P210325T3 <- ReadVisiumData(infoTable)
saveRDS(P210325T3, file = "P210325T3_HNSCC_semla.rds")

###Load the annotated scRNA-seq data for reference.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
SC_HNSCC <- readRDS(file = "HNSCC_Tumor_merge.rds") ##HNSCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of HNSCC tumors. The specific code is in the file named HNSCC_Tumor_merge.R.
table(SC_HNSCC$celltype)

###Normalize data and find variable features for ST data.
P210325T3 <- P210325T3 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000)

###Normalize scRNA-seq data and run vanilla analysis to create UMAP embedding.
SC_HNSCC <- SC_HNSCC |>
  NormalizeData() |>
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA() |> 
  RunUMAP(reduction = "pca", dims = 1:30)

###Rerun FindVariableFeatures to increase the number before cell type deconvolution.
SC_HNSCC <- SC_HNSCC |> 
  FindVariableFeatures(nfeatures = 10000)

###Visualize the available cell types in our UMAP embedding of the cells.
DimPlot(SC_HNSCC, group.by = "celltype")

###Run NNLS.(a quick method based on Non-Negative Least Squares to infer cell type proportions directly from Visium spot expression profiles.)
DefaultAssay(P210325T3) <- "Spatial"
ti <- Sys.time()
P210325T3 <- RunNNLS(object = P210325T3, 
                     singlecell_object = SC_HNSCC, 
                     groups = "celltype")
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))

###Check available cell types.
rownames(P210325T3)

###Selected cell types.
DefaultAssay(P210325T3) <- "celltypeprops"
selected_celltypes <- rownames(P210325T3)

###Loading H&E images.
setwd("E:/Raw Data/ST_data/HNSCC")
P210325T3 <- LoadImages(P210325T3, image_height = 1e3)
saveRDS(P210325T3, file = "P210325T3_HNSCC_semla2.rds")

###Visualization.
#Visualize single cell type.
setwd("E:/Raw Data/ST_data/HNSCC")
#MapFeatures plot for Supplementary Figure 10F: left.
pdf(file="ST_HNSCC_MapFeatures_P210325T3_Tumor_cells.pdf",width=4,height=5)
MapFeatures(P210325T3, pt_size = 2,
            features = selected_celltypes[13], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#MapFeatures plot for Supplementary Figure 10F: right.
pdf(file="ST_HNSCC_MapFeatures_P210325T3_apCAFs.pdf",width=4,height=5)
MapFeatures(P210325T3, pt_size = 2,
            features = selected_celltypes[14], image_use = "raw", section_number = 1,
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
setwd("E:/Raw Data/ST_data/HNSCC")
P210325T3 <- P210325T3 |> 
  LoadImages()

#MapMultipleFeatures plot for Supplementary Figure 10G.
pdf(file="ST_HNSCC_MapMultipleFeatures_P210325T3_Tumor_cells_apCAFs.pdf",width=6,height=4)
MapMultipleFeatures(P210325T3,
                    image_use = "raw",
                    section_number = 1,
                    label_by = "section_id",
                    ncol = 1,
                    pt_size = 2.1, max_cutoff = 0.99,
                    override_plot_dims = TRUE, 
                    colors = c("#88CCEE", "#CC6677"),
                    features = selected_celltypes[c(13, 14)])
dev.off()

###Obtain the selected cell type proportions for each spot in the tissue slice.
t <- FetchData(P210325T3, selected_celltypes[c(13, 14)])
names(t) <- c("Tumor_cells", "apCAFs")#rename cell type name.

###Correlation analysis.
data <- t
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor_cells)

###Visualization.
#scatter plot for Supplementary Figure 10H.
setwd("E:/Raw Data/ST_data/HNSCC")
pdf(file = "ST_HNSCC_P210325T3_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor_cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 10I, J.
###Load HNSCC Spatial transcriptome (ST) data.
#P210325T3 slice.
setwd("E:/Raw Data/ST_data/HNSCC/P210325T3")
t <- Read10X_h5("filtered_feature_bc_matrix.h5")
P210325T3 <- Seurat::CreateSeuratObject(counts = t, project = 'P210325T3', assay = 'Spatial')
P210325T3
P210325T3$slice <- "P210325T3"
S210325T3 <- Seurat::Read10X_Image(image.dir = './Spatial/', image.name = "tissue_hires_image.png")
S210325T3@scale.factors$lowres = S210325T3@scale.factors$hires
Seurat::DefaultAssay(object = S210325T3) <- 'Spatial'
S210325T3 <- S210325T3[colnames(x = P210325T3)]
P210325T3[['P210325T3']] <- S210325T3
P210325T3
setwd("E:/Raw Data/ST_data/HNSCC")
saveRDS(P210325T3, file = "HNSCC_ST_Seurat_P210325T3.rds")

###Apply sctransform normalization.
P210325T3  <- SCTransform(P210325T3 , assay = "Spatial", verbose = FALSE)
saveRDS(P210325T3, file = "HNSCC_ST_Seurat_P210325T3_sctransform.rds")

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/HNSCC")
gs = lapply(readLines("HNSCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
P210325T3 <- AddModuleScore(
  object = P210325T3,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
head(P210325T3@meta.data)
colnames(P210325T3@meta.data)[7] <- 'apCAFs signature'
colnames(P210325T3@meta.data)[8] <- 'CD4+ effector T cells signature'

##Visualization.
#SpatialFeaturePlot for Supplementary Figure 10I: left.
setwd("E:/Raw Data/ST_data/HNSCC")
pdf(file="HNSCC_ST_SpatialFeaturePlot_P210325T3_CD4_effector_T.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(P210325T3, features = "CD4+ effector T cells signature", pt.size.factor = 4.2, ncol = 1, crop = FALSE, images = "P210325T3")
dev.off()

#SpatialFeaturePlot for Supplementary Figure 10I: right.
pdf(file="HNSCC_ST_SpatialFeaturePlot_P210325T3_apCAFs.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(P210325T3, features = "apCAFs signature", pt.size.factor = 4.2, ncol = 1, crop = FALSE, images = "P210325T3")
dev.off()

###Obtain module scores for cell type signatures for each spot in the tissue slice.
m <- P210325T3@meta.data
df_P210325T3 <- m[, c(4, 7, 8)]

###Rename.
colnames(df_P210325T3) <- c("slice", "apCAFs", "CD4T")

###Correlation analysis.
###Shapiro-Wilk test for normality.
shapiro.test(df_P210325T3$CD4T)
shapiro.test(df_P210325T3$apCAFs)

###Visualization.
#scatter plot for Supplementary Figure 10J.
setwd("E:/Raw Data/ST_data/HNSCC")
pdf(file = "HNSCC_ST_P210325T3_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_P210325T3, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
