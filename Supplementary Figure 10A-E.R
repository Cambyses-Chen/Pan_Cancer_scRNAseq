###codes for Supplementary Figure 10A-E.

library(semla) 
library(tibble) 
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#1 Supplementary Figure 10A, B, C.
#HNSCC201125T10 slice.
###Load HNSCC Spatial transcriptome (ST) data.
#Read the file path.
setwd("E:/Raw Data/ST_data/HNSCC")
samples <- Sys.glob("./HNSCC201125T10/filtered_feature_bc_matrix.h5")
imgs <- Sys.glob("./HNSCC201125T10/tissue_hires_image.png")
spotfiles <- Sys.glob("./HNSCC201125T10/tissue_positions_list.csv")
json <- Sys.glob("./HNSCC201125T10/scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, 
                    section_id = c("HNSCC201125T10"))

###Create ST Seurat object.
HNSCC201125T10 <- ReadVisiumData(infoTable)
saveRDS(HNSCC201125T10, file = "HNSCC201125T10_HNSCC_semla.rds")

###Load the annotated scRNA-seq data for reference.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
SC_HNSCC <- readRDS(file = "HNSCC_Tumor_merge.rds") ##HNSCC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of HNSCC tumors. The specific code is in the file named HNSCC_Tumor_merge.R.
table(SC_HNSCC$celltype)

###Normalize data and find variable features for ST data.
HNSCC201125T10 <- HNSCC201125T10 |>
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
DefaultAssay(HNSCC201125T10) <- "Spatial"
ti <- Sys.time()
HNSCC201125T10 <- RunNNLS(object = HNSCC201125T10, 
                          singlecell_object = SC_HNSCC, 
                          groups = "celltype")
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))

###Check available cell types.
rownames(HNSCC201125T10)

###Selected cell types.
DefaultAssay(HNSCC201125T10) <- "celltypeprops"
selected_celltypes <- rownames(HNSCC201125T10)

###Loading H&E images.
setwd("E:/Raw Data/ST_data/HNSCC")
HNSCC201125T10 <- LoadImages(HNSCC201125T10, image_height = 1e3)
saveRDS(HNSCC201125T10, file = "HNSCC201125T10_HNSCC_semla2.rds")

###Visualization.
#Visualize single cell type.
setwd("E:/Raw Data/ST_data/HNSCC")
#MapFeatures plot for Supplementary Figure 10A: left.
pdf(file="ST_HNSCC_MapFeatures_HNSCC201125T10_Tumor_cells.pdf",width=4,height=5)
MapFeatures(HNSCC201125T10, pt_size = 2,
            features = selected_celltypes[13], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#MapFeatures plot for Supplementary Figure 10A: right.
pdf(file="ST_HNSCC_MapFeatures_HNSCC201125T10_apCAFs.pdf",width=4,height=5)
MapFeatures(HNSCC201125T10, pt_size = 2,
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
HNSCC201125T10 <- HNSCC201125T10 |> 
  LoadImages()

#MapMultipleFeatures plot for Supplementary Figure 10B.
pdf(file="ST_HNSCC_MapMultipleFeatures_HNSCC201125T10_Tumor_cells_apCAFs.pdf",width=6,height=4)
MapMultipleFeatures(HNSCC201125T10,
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
t <- FetchData(HNSCC201125T10, selected_celltypes[c(13, 14)])
names(t) <- c("Tumor_cells", "apCAFs")#rename cell type name.

###Correlation analysis.
data <- t
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor_cells)

###Visualization.
#scatter plot for Supplementary Figure 10C.
setwd("E:/Raw Data/ST_data/HNSCC")
pdf(file = "ST_HNSCC_HNSCC201125T10_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor_cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 10D, E.
###Load HNSCC Spatial transcriptome (ST) data.
#HNSCC201125T10 slice.
setwd("E:/Raw Data/ST_data/HNSCC/HNSCC201125T10")
t <- Read10X_h5("filtered_feature_bc_matrix.h5")
HNSCC201125T10 <- Seurat::CreateSeuratObject(counts = t, project = 'HNSCC201125T10', assay = 'Spatial')
HNSCC201125T10
HNSCC201125T10$slice <- "HNSCC201125T10"
S201125T10 <- Seurat::Read10X_Image(image.dir = './Spatial/', image.name = "tissue_hires_image.png")
S201125T10@scale.factors$lowres = S201125T10@scale.factors$hires
Seurat::DefaultAssay(object = S201125T10) <- 'Spatial'
S201125T10 <- S201125T10[colnames(x = HNSCC201125T10)]
HNSCC201125T10[['HNSCC201125T10']] <- S201125T10
HNSCC201125T10
setwd("E:/Raw Data/ST_data/HNSCC")
saveRDS(HNSCC201125T10, file = "HNSCC_ST_Seurat_HNSCC201125T10.rds")

###Apply sctransform normalization.
HNSCC201125T10  <- SCTransform(HNSCC201125T10 , assay = "Spatial", verbose = FALSE)
saveRDS(HNSCC201125T10, file = "HNSCC_ST_Seurat_HNSCC201125T10_sctransform.rds")

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/HNSCC")
gs = lapply(readLines("HNSCC_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
HNSCC201125T10 <- AddModuleScore(
  object = HNSCC201125T10,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
head(HNSCC201125T10@meta.data)
colnames(HNSCC201125T10@meta.data)[7] <- 'apCAFs signature'
colnames(HNSCC201125T10@meta.data)[8] <- 'CD4+ effector T cells signature'

###Visualization.
setwd("E:/Raw Data/ST_data/HNSCC")
#SpatialFeaturePlot for Supplementary Figure 10D: left.
pdf(file="HNSCC_ST_SpatialFeaturePlot_HNSCC201125T10_CD4_effector_T.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(HNSCC201125T10, features = "CD4+ effector T cells signature", pt.size.factor = 4.2, ncol = 1, crop = FALSE, images = "HNSCC201125T10")
dev.off()

#SpatialFeaturePlot for Supplementary Figure 10D: right.
pdf(file="HNSCC_ST_SpatialFeaturePlot_HNSCC201125T10_apCAFs.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(HNSCC201125T10, features = "apCAFs signature", pt.size.factor = 4.2, ncol = 1, crop = FALSE, images = "HNSCC201125T10")
dev.off()

###Obtain module scores for cell type signatures for each spot in the tissue slice.
m <- HNSCC201125T10@meta.data
df_HNSCC201125T10 <- m[, c(4, 7, 8)]

###Rename.
colnames(df_HNSCC201125T10) <- c("slice", "apCAFs", "CD4T")

###Correlation analysis.
#Shapiro-Wilk test for normality.
shapiro.test(df_HNSCC201125T10$CD4T)
shapiro.test(df_HNSCC201125T10$apCAFs)

###Visualization.
#scatter plot for Supplementary Figure 10E.
setwd("E:/Raw Data/ST_data/HNSCC")
pdf(file = "HNSCC_ST_HNSCC201125T10_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_HNSCC201125T10, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
