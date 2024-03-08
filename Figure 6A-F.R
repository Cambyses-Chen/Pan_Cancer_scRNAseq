###codes for Figure 6A-F.

library(Seurat)
library(scales)
library(ggsci)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(GSVA)
library(GSEABase)
library(dplyr)
library(msigdbr)

#1 Figure 6A.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

###Normalizing the data.
HNSCC_Fib <- NormalizeData(HNSCC_Fib, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
HNSCC_Fib <- FindVariableFeatures(HNSCC_Fib, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(HNSCC_Fib)
HNSCC_Fib <- ScaleData(HNSCC_Fib, features = all.genes)
HNSCC_Fib <- ScaleData(HNSCC_Fib, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
HNSCC_Fib <- RunPCA(HNSCC_Fib, features = VariableFeatures(object = HNSCC_Fib))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
HNSCC_Fib <- RunHarmony(HNSCC_Fib, group.by.vars = "orig.ident")
HNSCC_Fib <- RunUMAP(HNSCC_Fib, reduction = "harmony", dims = 1:20)
HNSCC_Fib <- RunTSNE(HNSCC_Fib, reduction = "harmony", dims = 1:20)
names(HNSCC_Fib@reductions)

###Quantify single-cell metabolism with Seurat.
HNSCC_Fib <- sc.metabolism.Seurat(obj = HNSCC_Fib, method = "AUCell", imputation = F, ncores = 6, metabolism.type = "KEGG")
saveRDS(HNSCC_Fib, file = "HNSCC_Fib_metabolism_KEGG.rds")

###Visualization.
#DotPlot for Figure 6A.
input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
pdf(file="HNSCC_Fib_DotPlot.metabolism_KEGG.pdf",width=6,height=5)
DotPlot.metabolism(obj = HNSCC_Fib, pathway = input.pathway, phenotype = "celltype", norm = "y") + xlab("") + theme(axis.text = element_text(size = 15)) + theme(axis.title.y = element_text(size = 15)) + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15))
dev.off()

#2 Figure 6B.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

###Normalizing the data.
NPC_Fib <- NormalizeData(NPC_Fib, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC_Fib <- FindVariableFeatures(NPC_Fib, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC_Fib)
NPC_Fib <- ScaleData(NPC_Fib, features = all.genes)
NPC_Fib <- ScaleData(NPC_Fib, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC_Fib <- RunPCA(NPC_Fib, features = VariableFeatures(object = NPC_Fib))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC_Fib <- RunHarmony(NPC_Fib, group.by.vars = "orig.ident")
NPC_Fib <- RunUMAP(NPC_Fib, reduction = "harmony", dims = 1:20)
NPC_Fib <- RunTSNE(NPC_Fib, reduction = "harmony", dims = 1:20)
names(NPC_Fib@reductions)

###Quantify single-cell metabolism with Seurat.
NPC_Fib <- sc.metabolism.Seurat(obj = NPC_Fib, method = "AUCell", imputation = F, ncores = 6, metabolism.type = "KEGG")
saveRDS(NPC_Fib, file = "NPC_Fib_metabolism_KEGG.rds")

###Visualization.
#DotPlot for Figure 6B.
input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
pdf(file="NPC_DotPlot.metabolism_KEGG.pdf",width=6,height=5)
DotPlot.metabolism(obj = NPC_Fib, pathway = input.pathway, phenotype = "celltype", norm = "y") + xlab("") + theme(axis.text = element_text(size = 15)) + theme(axis.title.y = element_text(size = 15)) + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15))
dev.off()

#3 Figure 6C.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_Fib$celltype)

###Normalizing the data.
BRCA_Fib <- NormalizeData(BRCA_Fib, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
BRCA_Fib <- FindVariableFeatures(BRCA_Fib, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(BRCA_Fib)
BRCA_Fib <- ScaleData(BRCA_Fib, features = all.genes)
BRCA_Fib <- ScaleData(BRCA_Fib, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
BRCA_Fib <- RunPCA(BRCA_Fib, features = VariableFeatures(object = BRCA_Fib))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
BRCA_Fib <- RunHarmony(BRCA_Fib, group.by.vars = "orig.ident")
BRCA_Fib <- RunUMAP(BRCA_Fib, reduction = "harmony", dims = 1:20)
BRCA_Fib <- RunTSNE(BRCA_Fib, reduction = "harmony", dims = 1:20)
names(BRCA_Fib@reductions)

###Quantify single-cell metabolism with Seurat.
BRCA_Fib <- sc.metabolism.Seurat(obj = BRCA_Fib, method = "AUCell", imputation = F, ncores = 6, metabolism.type = "KEGG")
saveRDS(BRCA_Fib, file = "BRCA_Fib_metabolism_KEGG.rds")

###Visualization.
#DotPlot for Figure 6C.
input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
pdf(file="BRCA_Fib_DotPlot.metabolism_KEGG.pdf",width=6,height=5)
DotPlot.metabolism(obj = BRCA_Fib, pathway = input.pathway, phenotype = "celltype", norm = "y") + xlab("") + theme(axis.text = element_text(size = 15)) + theme(axis.title.y = element_text(size = 15)) + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15))
dev.off()

#4 Figure 6D.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

###Extraction the normalization matrix.
expr <- as.matrix(HNSCC_Fib@assays$RNA@data)

###Load gene signature from msigdbr R package.
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
h_gene_sets_list = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

###Calculate GSVA scores for gene signatures in each cell.
gs.exp_h <- gsva(expr, h_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
write.csv(gs.exp_h, file = "gs.exp_HNSCC_Fib_H.csv")
saveRDS(gs.exp_h, file = "gs.exp_HNSCC_Fib_H.rds")

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
exp_g$barcode <- rownames(exp_g)

###Extract cell types information from the Seurat object.
ID <- HNSCC_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###HALLMARK_GLYCOLYSIS GSVA scores for gene signatures in each cell are integrated with cell types information.
exp_g2 <- merge(ID2, exp_g, by.x = "barcode", by.y = "barcode")
rownames(exp_g2) <- exp_g2[, 1]

###Anderson-Darling test for normality.
library(nortest)
ad.test(exp_g2$HALLMARK_GLYCOLYSIS)

###Visualization.
library(ggpubr)
table(exp_g2$celltype)

#boxplot for Figure 6D.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"), c("apCAFs", "NFs") )
pdf(file = "ggboxplot_HNSCC_GLYCOLYSIS.pdf", width=5, height=4.5)
ggboxplot(exp_g2, x = "celltype", y = "HALLMARK_GLYCOLYSIS", color = "celltype", 
          add = "jitter", legend = "none", title = "HALLMARK GLYCOLYSIS", xlab = "", ylab = "GSVA enrichment score", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.2, 0.3, 0.4, 0.5))
dev.off()

#5 Figure 6E.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

###Extraction the normalization matrix.
expr <- as.matrix(NPC_Fib@assays$RNA@data)

###Load gene signature from msigdbr R package.
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
h_gene_sets_list = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

###Calculate GSVA scores for gene signatures in each cell.
gs.exp_h <- gsva(expr, h_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
write.csv(gs.exp_h, file = "gs.exp_NPC_Fib_H.csv")
saveRDS(gs.exp_h, file = "gs.exp_NPC_Fib_H.rds")

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
exp_g$barcode <- rownames(exp_g)

###Extract cell types information from the Seurat object.
ID <- NPC_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###HALLMARK_GLYCOLYSIS GSVA scores for gene signatures in each cell are integrated with cell types information.
exp_g2 <- merge(ID2, exp_g, by.x = "barcode", by.y = "barcode")
rownames(exp_g2) <- exp_g2[, 1]

###Anderson-Darling test for normality.
library(nortest)
ad.test(exp_g2$HALLMARK_GLYCOLYSIS)

###Visualization.
library(ggpubr)
table(exp_g2$celltype)

#boxplot for Figure 6E.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"), c("apCAFs", "NFs") )
pdf(file = "ggboxplot_NPC_GLYCOLYSIS.pdf", width=5, height=4.5)
ggboxplot(exp_g2, x = "celltype", y = "HALLMARK_GLYCOLYSIS", color = "celltype", 
          add = "jitter", legend = "none", title = "HALLMARK GLYCOLYSIS", xlab = "", ylab = "GSVA enrichment score", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.2, 0.3, 0.4, 0.5))
dev.off()

#6 Figure 6F.
#BRCA.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BRCA_Fib$celltype)

###Extraction the normalization matrix.
expr <- as.matrix(BRCA_Fib@assays$RNA@data)

###Load gene signature from msigdbr R package.
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
h_gene_sets_list = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

###Calculate GSVA scores for gene signatures in each cell.
gs.exp_h <- gsva(expr, h_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
write.csv(gs.exp_h, file = "gs.exp_BRCA_Fib_H.csv")
saveRDS(gs.exp_h, file = "gs.exp_BRCA_Fib_H.rds")

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
exp_g$barcode <- rownames(exp_g)

###Extract cell types information from the Seurat object.
ID <- BRCA_Fib@meta.data
ID$barcode <- rownames(ID)
ID2 <- ID[, c(11, 10)]

###HALLMARK_GLYCOLYSIS GSVA scores for gene signatures in each cell are integrated with cell types information.
exp_g2 <- merge(ID2, exp_g, by.x = "barcode", by.y = "barcode")
rownames(exp_g2) <- exp_g2[, 1]

###Anderson-Darling test for normality.
library(nortest)
ad.test(exp_g2$HALLMARK_GLYCOLYSIS)

###Visualization.
library(ggpubr)
table(exp_g2$celltype)

#boxplot for Figure 6F.
my_comparisons <- list( c("apCAFs", "myCAFs"), c("apCAFs", "eCAFs"), c("apCAFs", "iCAFs"))
pdf(file = "ggboxplot_BRCA_GLYCOLYSIS.pdf", width=5, height=4.5)
ggboxplot(exp_g2, x = "celltype", y = "HALLMARK_GLYCOLYSIS", color = "celltype", 
          add = "jitter", legend = "none", title = "HALLMARK GLYCOLYSIS", xlab = "", ylab = "GSVA enrichment score", font.label = list(size = 20, face = "bold")) +
  rotate_x_text(angle = 45) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", label.y = c(0.2, 0.3, 0.4, 0.5))
dev.off()
