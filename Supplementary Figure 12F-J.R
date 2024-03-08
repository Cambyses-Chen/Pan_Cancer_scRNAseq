###codes for Supplementary Figure 12F-J.

library(semla) 
library(tibble) 
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(Matrix)

#1 Supplementary Figure 12F, G, H.
#CID4535 slice.
###Load BRCA Spatial transcriptome (ST) data.
#Read the file path.
setwd("E:/Raw Data/ST_data/BRCA")
samples <- Sys.glob("./CID4535_spatial/filtered_count_matrices_folder")
imgs <- Sys.glob("./CID4535_spatial/Spatial/tissue_hires_image.png")
spotfiles <- Sys.glob("./CID4535_spatial/Spatial/tissue_positions_list.csv")
json <- Sys.glob("./CID4535_spatial/Spatial/scalefactors_json.json")
infoTable <- tibble(samples, imgs, spotfiles, json, 
                    section_id = c("CID4535"))

###Create ST Seurat object.
CID4535 <- ReadVisiumData(infoTable)
saveRDS(CID4535, file = "CID4535_BRCA_semla.rds")

###Load the annotated scRNA-seq data for reference.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
SC_BC <- readRDS(file = "BC_Tumor_merge.rds") #BC_Tumor_merge.rds is a Seurat object annotated with various cell types for single-cell RNA-seq samples of BRCA tumors. The specific code is in the file named BC_Tumor_merge.R.
table(SC_BC$celltype)

###Normalize data and find variable features for ST data.
CID4535 <- CID4535 |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000)

###Normalize scRNA-seq data and run vanilla analysis to create UMAP embedding.
SC_BC <- SC_BC |>
  NormalizeData() |>
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA() |> 
  RunUMAP(reduction = "pca", dims = 1:30)

###Rerun FindVariableFeatures to increase the number before cell type deconvolution.
SC_BC <- SC_BC |> 
  FindVariableFeatures(nfeatures = 10000)

###Visualize the available cell types in our UMAP embedding of the cells.
DimPlot(SC_BC, group.by = "celltype")

###Run NNLS.(a quick method based on Non-Negative Least Squares to infer cell type proportions directly from Visium spot expression profiles.)
DefaultAssay(CID4535) <- "Spatial"
ti <- Sys.time()
CID4535 <- RunNNLS(object = CID4535, 
                   singlecell_object = SC_BC, 
                   groups = "celltype")
sprintf("RunNNLS completed in %s seconds", round(Sys.time() - ti, digits = 2))

###Check available cell types.
rownames(CID4535)

###Selected cell types.
DefaultAssay(CID4535) <- "celltypeprops"
selected_celltypes <- rownames(CID4535)

###Loading H&E images.
setwd("E:/Raw Data/ST_data/BRCA")
CID4535 <- LoadImages(CID4535, image_height = 1e3)
saveRDS(CID4535, file = "CID4535_BRCA_semla2.rds")

###Visualization.
#Visualize single cell type.
setwd("E:/Raw Data/ST_data/BRCA")
#MapFeatures plot for Supplementary Figure 12F: left.
pdf(file="ST_BRCA_MapFeatures_CID4535_Tumor_cells.pdf",width=4,height=5)
MapFeatures(CID4535, pt_size = 1.65,
            features = selected_celltypes[11], image_use = "raw", section_number = 1,
            arrange_features = "row", scale = "shared", label_by = "section_id",
            override_plot_dims = TRUE, ncol = 1,
            colors = RColorBrewer::brewer.pal(n = 9, name = "Spectral") |> rev(), 
            scale_alpha = TRUE)  + 
  theme(legend.position = "right", legend.margin = margin(b = 50),
        legend.text = element_text(angle = 0),
        plot.title = element_blank())
dev.off()

#MapFeatures plot for Supplementary Figure 12F: right.
pdf(file="ST_BRCA_MapFeatures_CID4535_apCAFs.pdf",width=4,height=5)
MapFeatures(CID4535, pt_size = 1.65,
            features = selected_celltypes[13], image_use = "raw", section_number = 1,
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
setwd("E:/Raw Data/ST_data/BRCA")
CID4535 <- CID4535 |> 
  LoadImages()

#MapMultipleFeatures plot for Supplementary Figure 12G.
pdf(file="ST_BRCA_MapMultipleFeatures_CID4535_Tumor_cells_apCAFs.pdf",width=6,height=4)
MapMultipleFeatures(CID4535,
                    image_use = "raw",
                    section_number = 1,
                    label_by = "section_id",
                    ncol = 1,
                    pt_size = 1.8, max_cutoff = 0.99,
                    override_plot_dims = TRUE, 
                    colors = c("#88CCEE", "#CC6677"),
                    features = selected_celltypes[c(11, 13)])
dev.off()

###Obtain the selected cell type proportions for each spot in the tissue slice.
t <- FetchData(CID4535, selected_celltypes[c(11, 13)])
names(t) <- c("Tumor_cells", "apCAFs")#rename cell type name.

###Correlation analysis.
data <- t
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$Tumor_cells)

###Visualization.
#Scatter plot for Supplementary Figure 12H.
setwd("E:/Raw Data/ST_data/BRCA")
pdf(file = "ST_BRCA_CID4535_apCAFs_Tumor_spearman.pdf", width=6, height=5)
ggscatter(data, x = "Tumor_cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tumor cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 12I, J.
###Load BRCA Spatial transcriptome (ST) data.
#CID4535 slice.
setwd("E:/Raw Data/ST_data/BRCA/CID4535_spatial")
matrix_dir = "./filtered_count_matrices_folder/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
t <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(t) = barcode.names$V1
rownames(t) = feature.names$V1
CID4535 <- Seurat::CreateSeuratObject(counts = t, project = 'CID4535', assay = 'Spatial')
CID4535
CID4535$slice <- "CID4535"
S4535 <- Seurat::Read10X_Image(image.dir = './Spatial/')
Seurat::DefaultAssay(object = S4535) <- 'Spatial'
S4535 <- S4535[colnames(x = CID4535)]
CID4535[['CID4535']] <- S4535
CID4535
setwd("E:/Raw Data/ST_data/BRCA")
saveRDS(CID4535, file = "BRCA_ST_Seurat_CID4535.rds")

###Apply sctransform normalization.
CID4535  <- SCTransform(CID4535 , assay = "Spatial", verbose = FALSE)
saveRDS(CID4535, file = "BRCA_ST_Seurat_CID4535_sctransform.rds")

###Load gene signature. apCAFs represents the apCAFs signature, CD4_T represents the CD4+ effector T cells signature.
setwd("E:/Raw Data/Signatures/BRCA")
gs = lapply(readLines("BRCA_apCAFs_CD4.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###AddModuleScore.
CID4535 <- AddModuleScore(
  object = CID4535,
  features = gs,
  ctrl = 100,
  name = c("apCAFs", "CD4T")
)

###Rename.
head(CID4535@meta.data)
colnames(CID4535@meta.data)[7] <- 'apCAFs signature'
colnames(CID4535@meta.data)[8] <- 'CD4+ effector T cells signature'

##Visualization.
setwd("E:/Raw Data/ST_data/BRCA")
#SpatialFeaturePlot for Supplementary Figure 12I: left.
pdf(file="BRCA_ST_SpatialFeaturePlot_CID4535_CD4_effector_T.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(CID4535, features = "CD4+ effector T cells signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "CID4535")
dev.off()

#SpatialFeaturePlot for Supplementary Figure 12I: right.
pdf(file="BRCA_ST_SpatialFeaturePlot_CID4535_apCAFs.pdf",width=3.5,height=4.2)
SpatialFeaturePlot(CID4535, features = "apCAFs signature", pt.size.factor = 1.3, ncol = 1, crop = FALSE, images = "CID4535")
dev.off()

###Obtain module scores for cell type signatures for each spot in the tissue slice.
m <- CID4535@meta.data
df_CID4535 <- m[, c(4, 7, 8)]

###Rename.
colnames(df_CID4535) <- c("slice", "apCAFs", "CD4T")

###Correlation analysis.
###Shapiro-Wilk test for normality.
shapiro.test(df_CID4535$CD4T)
shapiro.test(df_CID4535$apCAFs)

###Visualization.
setwd("E:/Raw Data/ST_data/BRCA")
#scatter plot for Supplementary Figure 12J.
pdf(file = "BRCA_ST_CID4535_apCAFs_CD4_effector_T_cells_spearman.pdf", width=6, height=5)
ggscatter(df_CID4535, x = "CD4T", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD4+ effector T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
