###codes for Figure 5A-F.

library(Seurat)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(ComplexHeatmap)
library(pheatmap)

#1 Figure 5A.
#HNSCC.
###R
###Extraction expression matrix.
###Load single-cell RNA sequencing data of Fibroblasts major cell type in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds")  #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/HNSCC")
write.csv(t(as.matrix(HNSCC_Fib@assays$RNA@counts)),file = "HNSCC_Fib_0408.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/HNSCC")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("HNSCC_Fib_0408.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("HNSCC_Fib_0408.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output HNSCC_Fib_0408_adj.sample.tsv \
--method grnboost2 \
HNSCC_Fib_0408.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
HNSCC_Fib_0408_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname HNSCC_Fib_0408.loom \
--mode "dask_multiprocessing" \
--output HNSCC_Fib_0408_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
HNSCC_Fib_0408.loom \
HNSCC_Fib_0408_reg.csv \
--output HNSCC_Fib_0408_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/HNSCC")
scenicLoomPath='HNSCC_Fib_0408_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
###Load single-cell RNA sequencing data of Fibroblasts major cell type in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib <- readRDS(file = "HNSCC_Fib_new.rds")  #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(HNSCC_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
HNSCC_Fib2 <- subset(HNSCC_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(HNSCC_Fib2))
cellInfo <- HNSCC_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(4, 3, 1, 2, 5)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/HNSCC")
write.csv(rss_data, file = "rss_data_HNSCC_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(HNSCC_Fib2), 
                        celltype = HNSCC_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("HNSCC_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
hnscc <- gs[[1]]
hnscc

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[hnscc, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("iCAFs",1),
                               rep("NFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC")
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "iCAFs", "NFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5A.
pdf(file="pheatmap_HNSCC_Fib_AUC.pdf",width=5, height=6)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "HNSCC", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.725,y=0.945, gp=gpar(fontsize=10, fontface="bold"))
dev.off()

#2 Figure 5B.
#OV.
###R
###Extraction expression matrix.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in OV.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV_Fib <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The codes for OV_Fib_res0.2_annotation_sub.rds is located at line 406 of Figure 1G-J.R.
table(OV_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/OV")
write.csv(t(as.matrix(OV_Fib@assays$RNA@counts)),file = "OV_Fib_0525.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/OV")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("OV_Fib_0525.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("OV_Fib_0525.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output OV_Fib_0525_adj.sample.tsv \
--method grnboost2 \
OV_Fib_0525.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
OV_Fib_0525_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname OV_Fib_0525.loom \
--mode "dask_multiprocessing" \
--output OV_Fib_0525_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
OV_Fib_0525.loom \
OV_Fib_0525_reg.csv \
--output OV_Fib_0525_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/OV")
scenicLoomPath='OV_Fib_0525_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in OV.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV_Fib <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The codes for OV_Fib_res0.2_annotation_sub.rds is located at line 406 of Figure 1G-J.R.
table(OV_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(OV_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
OV_Fib2 <- subset(OV_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(OV_Fib2))
cellInfo <- OV_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(4, 2, 1, 3)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/OV")
write.csv(rss_data, file = "rss_data_OV_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(OV_Fib2), 
                        celltype = OV_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("OV_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
ov <- gs[[1]]
ov

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[ov, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("iCAFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574")
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "iCAFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5B.
pdf(file="pheatmap_OV_Fib_AUC.pdf",width=5, height=5)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "OV", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.69,y=0.93, gp=gpar(fontsize=10, fontface="bold"))
dev.off()

#3 Figure 5C.
#NPC.
###R
###Extraction expression matrix.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/NPC")
write.csv(t(as.matrix(NPC_Fib@assays$RNA@counts)),file = "NPC_Fib_0408.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/NPC")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("NPC_Fib_0408.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("NPC_Fib_0408.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output NPC_Fib_0408_adj.sample.tsv \
--method grnboost2 \
NPC_Fib_0408.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
NPC_Fib_0408_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname NPC_Fib_0408.loom \
--mode "dask_multiprocessing" \
--output NPC_Fib_0408_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
NPC_Fib_0408.loom \
NPC_Fib_0408_reg.csv \
--output NPC_Fib_0408_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/NPC")
scenicLoomPath='NPC_Fib_0408_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(NPC_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
NPC_Fib2 <- subset(NPC_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(NPC_Fib2))
cellInfo <- NPC_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(1, 3, 2, 4, 5)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/NPC")
write.csv(rss_data, file = "rss_data_NPC_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(NPC_Fib2), 
                        celltype = NPC_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("NPC_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
npc <- gs[[1]]
npc

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[npc, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("iCAFs",1),
                               rep("NFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC")
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "iCAFs", "NFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5C.
pdf(file="pheatmap_NPC_Fib_AUC.pdf",width=5, height=6)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "NPC", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.719,y=0.945, gp=gpar(fontsize=10, fontface="bold"))
dev.off()

#4 Figure 5D.
#BRCA.
###R
###Extraction expression matrix.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BC_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BC_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/BRCA")
write.csv(t(as.matrix(BC_Fib@assays$RNA@counts)),file = "BC_Fib_0408.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/BRCA")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("BC_Fib_0408.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("BC_Fib_0408.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output BC_Fib_0408_adj.sample.tsv \
--method grnboost2 \
BC_Fib_0408.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
BC_Fib_0408_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname BC_Fib_0408.loom \
--mode "dask_multiprocessing" \
--output BC_Fib_0408_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
BC_Fib_0408.loom \
BC_Fib_0408_reg.csv \
--output BC_Fib_0408_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/BRCA")
scenicLoomPath='BC_Fib_0408_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in BRCA.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BC_Fib <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The codes for BC_GSE176078_Fib_res0.2_annotation2.rds is located at line 437 of Supplementary Figure 5G-H.R.
table(BC_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(BC_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
BC_Fib2 <- subset(BC_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(BC_Fib2))
cellInfo <- BC_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(2, 3, 4, 1)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/BRCA")
write.csv(rss_data, file = "rss_data_BRCA_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(BC_Fib2), 
                        celltype = BC_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("BRCA_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
brca <- gs[[1]]
brca

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[brca, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("iCAFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574")
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "iCAFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5D.
pdf(file="pheatmap_BRCA_Fib_AUC.pdf",width=5, height=5)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "BRCA", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.69,y=0.93, gp=gpar(fontsize=10, fontface="bold"))
dev.off()

#5 Figure 5E.
#CM.
###R
###Extraction expression matrix.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in CM.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM_Fib <- readRDS(file = "CM_Fib_res0.4_annotation_sub.rds") #The codes for CM_Fib_res0.4_annotation_sub.rds is located at line 392 of Supplementary Figure 5K-L.R.
table(CM_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/CM")
write.csv(t(as.matrix(CM_Fib@assays$RNA@counts)),file = "CM_Fib_0408.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/CM")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("CM_Fib_0408.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("CM_Fib_0408.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output CM_Fib_0408_adj.sample.tsv \
--method grnboost2 \
CM_Fib_0408.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
CM_Fib_0408_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname CM_Fib_0408.loom \
--mode "dask_multiprocessing" \
--output CM_Fib_0408_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
CM_Fib_0408.loom \
CM_Fib_0408_reg.csv \
--output CM_Fib_0408_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/CM")
scenicLoomPath='CM_Fib_0408_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
#Load annotated single-cell RNA sequencing data of CAFs subpopulations in CM.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM_Fib <- readRDS(file = "CM_Fib_res0.4_annotation_sub.rds") #The codes for CM_Fib_res0.4_annotation_sub.rds is located at line 392 of Supplementary Figure 5K-L.R.
table(CM_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(CM_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
CM_Fib2 <- subset(CM_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(CM_Fib2))
cellInfo <- CM_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(3, 2, 4, 1)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/CM")
write.csv(rss_data, file = "rss_data_CM_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(CM_Fib2), 
                        celltype = CM_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("CM_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
cm <- gs[[1]]
cm

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[cm, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("iCAFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574")
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "iCAFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5E.
pdf(file="pheatmap_CM_Fib_AUC.pdf",width=5, height=5)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "CM", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.685,y=0.93, gp=gpar(fontsize=10, fontface="bold"))
dev.off()

#6 Figure 5F.
#RCC.
###R
###Extraction expression matrix.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in RCC.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC_Fib <- readRDS(file = "RCC_Fib_new.rds") #The codes for RCC_Fib_new.rds is located at line 159 of Supplementary Figure 9F-G.R.
table(RCC_Fib$celltype)

setwd("E:/Raw Data/pySCENIC_data/RCC")
write.csv(t(as.matrix(RCC_Fib@assays$RNA@counts)),file = "RCC_Fib_0525.csv")

###python3
###Convert csv file to loom file.
import os
os.chdir("E:/Raw Data/pySCENIC_data/RCC")
os.getcwd()
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("RCC_Fib_0525.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("RCC_Fib_0525.loom",x.X.transpose(),row_attrs,col_attrs);

###linux pyscenic conda environment.
###pySCENIC analysis.
pyscenic grn \
--num_workers 10 \
--output RCC_Fib_0525_adj.sample.tsv \
--method grnboost2 \
RCC_Fib_0525.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
RCC_Fib_0525_adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname RCC_Fib_0525.loom \
--mode "dask_multiprocessing" \
--output RCC_Fib_0525_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
RCC_Fib_0525.loom \
RCC_Fib_0525_reg.csv \
--output RCC_Fib_0525_SCENIC.loom \
--num_workers 10

###R
###Read the output file.
setwd("E:/Raw Data/pySCENIC_data/RCC")
scenicLoomPath='RCC_Fib_0525_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

###Get AUCell matrix from the loom file.
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
dim(regulonAUC)

###Extract cell types information from the Seurat object.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in RCC.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC_Fib <- readRDS(file = "RCC_Fib_new.rds") #The codes for RCC_Fib_new.rds is located at line 159 of Supplementary Figure 9F-G.R.
table(RCC_Fib$celltype)

###The AUCell submatrix corresponding to the cell barcodes was acquired.
common <- intersect(colnames(RCC_Fib), colnames(regulonAUC))
sub_regulonAUC <- regulonAUC[, common]
dim(sub_regulonAUC)
RCC_Fib2 <- subset(RCC_Fib, cells = common)
identical(colnames(sub_regulonAUC), colnames(RCC_Fib2))
cellInfo <- RCC_Fib2@meta.data
table(cellInfo$celltype)

###Computes the regulon specificity score(rss) for each cell type.
rss <- calcRSS(AUC=getAUC(sub_regulonAUC), cellAnnotation=cellInfo$celltype)
rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 0,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
rss_data <- rssPlot$plot$data

###According to the rss, we got Z-scores of regulon for each cell type.
library(reshape2)
rss_data <- dcast(rss_data, 
                  Topic~rss_data$cellType,
                  value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
rss_data <- rss_data[,c(1, 3, 2, 4)]
colnames(rss_data)
setwd("E:/Raw Data/pySCENIC_data/RCC")
write.csv(rss_data, file = "rss_data_RCC_Z.csv")

###Select the regulons that are enriched specifically for individual cell types.

###Generate a matrix of AUC scores corresponding to regulons for every cell type.
cellTypes <- data.frame(row.names = colnames(RCC_Fib2), 
                        celltype = RCC_Fib2$celltype)
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes$celltype)
regulonActivity <- sapply(cellsPerGroup,
                          function(cells)rowMeans(getAUC(sub_regulonAUC)[,cells]))###Calculate average expression.

###According to the AUC scores, we got Z-scores of regulon for each cell type.
regulonActivity_Scaled <- t(scale(t(regulonActivity), center = T, scale=T))
dim(regulonActivity_Scaled)
regulonActivity_Scaled=na.omit(regulonActivity_Scaled)
colnames(regulonActivity_Scaled)

###Read the regulons that are enriched specifically for individual cell types.
gs = lapply(readLines("RCC_regulons.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs
rcc <- gs[[1]]
rcc

###AUC z-score submatrix of regulons who have been chosen.
regulonActivity_Scaled <- regulonActivity_Scaled[rcc, ]

###Visualization.
col_ann <- data.frame(group= c(rep("apCAFs",1),
                               rep("myCAFs",1),
                               rep("eCAFs",1),
                               rep("NFs",1)))
rownames(col_ann) <- colnames(regulonActivity_Scaled)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#0099CC") 
names(groupcol) <- c("apCAFs", "myCAFs", "eCAFs", "NFs")
col <- list(group=groupcol)

#heatmap plot for Figure 5F.
pdf(file="pheatmap_RCC_Fib_AUC.pdf",width=5, height=5)
pheatmap(regulonActivity_Scaled, color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_rows = T,cluster_cols = F,
         annotation_col = col_ann,
         annotation_color = col, legend = TRUE,
         show_colnames = FALSE,
         fontsize = 10, main = "RCC", cellwidth = 25, cellheight = 25)
grid.text("Z-score", x=0.68,y=0.93, gp=gpar(fontsize=10, fontface="bold"))
dev.off()
