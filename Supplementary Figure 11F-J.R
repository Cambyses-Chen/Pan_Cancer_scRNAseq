###codes for Supplementary Figure 11F-J.

library(semla) 
library(tibble) 
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#1 Supplementary Figure 11F, G, H.
#GSM6592137 slice.
###Load OV Spatial transcriptome (ST) data.
#Read the file path.
setwd("E:/Raw Data/ST_data/OV")
samples <- Sys.glob("./GSM6592137/filtered_count_matrices_folder")
imgs <- Sys.glob("./GSM6592137/tissue_hires_image.png")
spotfiles <- Sys.glob("./GSM6592137/tissue_positions_list.csv")
json <- Sys.glob("./GSM6592137/scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, 
                    section_id = c("GSM6592137"))

###Create ST Seurat object.
GSM6592137 <- ReadVisiumData(infoTable)
saveRDS(GSM6592137, file = "GSM6592137_OV_semla.rds")

###Load the annotated scRNA-seq data for reference.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
SC_OV <- readRDS(file = "OV_naive_merge.rds") ##OV_naive_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of OV tumors without drug treatment. The specific code is in the file named OV_naive_merge.R.
table(SC_OV$celltype)

###Normalize data and find variable features for ST data.
GSM6592137 <- GSM6592137 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000)

###Normalize scRNA-seq data and run vanilla analysis to create UMAP embedding.
SC_OV <- SC_OV |>
  NormalizeData() |>
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA() |> 
  RunUMAP(reduction = "pca", dims = 1:30)

###Rerun FindVariableFeatures to increase the number before cell type deconvolution.
SC_OV <- SC_OV |> 
  FindVariableFeatures(nfeatures = 10000)

###Visualize the available cell types in our UMAP embedding of the cells.
DimPlot(SC_OV, group.by = "celltype")

###Run NNLS.(a quick method based on Non-Negative Least Squares to infer cell type proportions directly from Visium spot expression profiles.)
DefaultAssay(GSM6592137) <- "Spatial"
ti <- Sys.time()
GSM6592137 <- RunNNLS(object = GSM6592137, 
                      singlecell_object = SC_OV, 
                      groups = "celltype")
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))

###Check available cell types.
rownames(GSM6592137)

###Selected cell types.
DefaultAssay(GSM6592137) <- "celltypeprops"
selected_celltypes <- rownames(GSM6592137)

###Loading H&E images.
setwd("E:/Raw Data/ST_data/OV")
GSM6592137 <- LoadImages(GSM6592137, image_height = 1e3)
saveRDS(GSM6592137, file = "GSM6592137_OV_semla2.rds")

###Visualization.
#Visualize single cell type.
setwd("E:/Raw Data/ST_data/OV")
#MapFeatures plot for Supplementary Figure 11F: left.
pdf(file="ST_OV_MapFeatures_GSM6592137_Tumor_cells.pdf",width=4,height=5)
MapFeatures(GSM6592137, pt_size = 1.7,
            features = selected_celltypes[9], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#MapFeatures plot for Supplementary Figure 11F: right.
pdf(file="ST_OV_MapFeatures_GSM6592137_apCAFs.pdf",width=4,height=5)
MapFeatures(GSM6592137, pt_size = 1.7,
            features = selected_celltypes[11], image_use = "raw", section_number = 1,
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
setwd("E:/Raw Data/ST_data/OV")
GSM6592137 <- GSM6592137 |> 
  LoadImages()

#MapMultipleFeatures plot for Supplementary Figure 11G.
pdf(file="ST_OV_MapMultipleFeatures_GSM6592137_Tumor_cells_apCAFs.pdf",width=6,height=4)
MapMultipleFeatures(GSM6592137,
                    image_use = "raw",
                    section_number = 1,
                    label_by = "section_id",
                    ncol = 1,
                    pt_size = 1.8, max_cutoff = 0.99,
                    override_plot_dims = TRUE, 
                    colors = c("#88CCEE", "#CC6677"),
                    features = selected_celltypes[c(9, 11)])
dev.off()

###Obtain the selected cell type proportions for each spot in the tissue slice.
t <- FetchData(GSM6592137, selected_celltypes[c(9, 11)])
names(t) <- c("Tumor_cells", "apCAFs")#rename cell type name.

###Correlation analysis.
data <- t
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor_cells)

###Visualization.
#Scatter plot for Supplementary Figure 11H.
setwd("E:/Raw Data/ST_data/OV")
pdf(file = "ST_OV_GSM6592137_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor_cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 11I, J.
###Load OV Spatial transcriptome (ST) data.
#GSM6592137 slice.
setwd("E:/Raw Data/ST_data/OV/GSM6592137")
t <- read.csv(file = "raw_counts.csv")
rownames(t) <- t[, 1]
t <- t[, -1]
names(t) <- gsub(x = names(t), pattern = "\\.", replacement = "-")
GSM6592137 <- Seurat::CreateSeuratObject(counts = t, project = 'GSM6592137', assay = 'Spatial')
GSM6592137$slice <- "GSM6592137"
S6592137 <- Seurat::Read10X_Image(image.dir = './Spatial/')
Seurat::DefaultAssay(object = S6592137) <- 'Spatial'
S6592137 <- S6592137[colnames(x = GSM6592137)]
GSM6592137[['GSM6592137']] <- S6592137
GSM6592137
setwd("E:/Raw Data/ST_data/OV")
saveRDS(GSM6592137, file = "OV_ST_Seurat_GSM6592137.rds")

###Apply sctransform normalization.
GSM6592137  <- SCTransform(GSM6592137 , assay = "Spatial", verbose = FALSE)
saveRDS(GSM6592137, file = "OV_ST_Seurat_GSM6592137_sctransform.rds")

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/OV")
gs = lapply(readLines("OV_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
GSM6592137 <- AddModuleScore(
  object = GSM6592137,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
head(GSM6592137@meta.data)
colnames(GSM6592137@meta.data)[7] <- 'apCAFs signature'
colnames(GSM6592137@meta.data)[8] <- 'CD4+ effector T cells signature'

###Visualization.
setwd("E:/Raw Data/ST_data/OV")
#SpatialFeaturePlot for Supplementary Figure 11I: left.
pdf(file="OV_ST_SpatialFeaturePlot_GSM6592137_CD4_effector_T.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(GSM6592137, features = "CD4+ effector T cells signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "GSM6592137")
dev.off()

#SpatialFeaturePlot for Supplementary Figure 11I: right.
pdf(file="OV_ST_SpatialFeaturePlot_GSM6592137_apCAFs.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(GSM6592137, features = "apCAFs signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "GSM6592137")
dev.off()

###Obtain module scores for cell type signatures for each spot in the tissue slice.
m <- GSM6592137@meta.data
df_GSM6592137 <- m[, c(4, 7, 8)]

###Rename.
colnames(df_GSM6592137) <- c("slice", "apCAFs", "CD4T")

###Correlation analysis.
###Shapiro-Wilk test for normality.
shapiro.test(df_GSM6592137$CD4T)
shapiro.test(df_GSM6592137$apCAFs)

###Visualization.
#scatter plot for Supplementary Figure 11J.
setwd("E:/Raw Data/ST_data/OV")
pdf(file = "OV_ST_GSM6592137_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_GSM6592137, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
