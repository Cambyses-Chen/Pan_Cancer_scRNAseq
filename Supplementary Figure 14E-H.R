###codes for Supplementary Figure 14E-H.

library(Seurat)
library(scales)
library(ggplot2)
library(corrplot)
library(ggsci)
library(reshape2)
library(monocle)
library(ggridges)

#1 Supplementary Figure 14E.
#RCC.
###Load Seurat objects of various annotated cell types from RCC single-cell RNA sequencing data.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC <- readRDS(file = "RCC_resolution0.1_annotation2.rds") #The codes for RCC_resolution0.1_annotation2.rds is located at line 288 of Supplementary Figure 5M-N.R.
RCC_M <- readRDS(file = "RCC_Myeloid_res0.1_annotation2.rds") #The codes for RCC_Myeloid_res0.1_annotation2.rds is located at line 176 of Supplementary Figure 6Q-R.R.
RCC_NF_ap <- readRDS(file = "RCC_NF_ap.rds") #The codes for RCC_NF_ap.rds is located at line 172 of Supplementary Figure 9F-G.R.
table(RCC$Celltype)
table(RCC_M$celltype)
table(RCC_NF_ap$celltype)

###Extract endothelial cells.
RCC_Endo <- subset(RCC, subset=Celltype %in% c("Endothelial cells"))
RCC_Endo$Celltype <- as.factor(as.character(RCC_Endo$Celltype))
table(RCC_Endo$Celltype)
RCC_Endo$celltype <- RCC_Endo$Celltype

###Extract the required myeloid cells.
RCC_M2 <- subset(RCC_M, subset=celltype %in% c("cDC1", "cDC2", "C1QC+ Macrophages", "SPP1+ Macrophages", "VCAN+ Macrophages"))
RCC_M2$celltype <- as.factor(as.character(RCC_M2$celltype))
table(RCC_M2$celltype)

###Merge into a Seurat object containing various antigen-presenting cells (APCs).
RCC_APC <- merge(x = RCC_Endo, y = c(RCC_M2, RCC_NF_ap))
table(RCC_APC$celltype)

###Normalizing the data.
RCC_APC <- NormalizeData(RCC_APC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
RCC_APC <- FindVariableFeatures(RCC_APC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(RCC_APC)
RCC_APC <- ScaleData(RCC_APC, features = all.genes)
RCC_APC <- ScaleData(RCC_APC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
RCC_APC <- RunPCA(RCC_APC, features = VariableFeatures(object = RCC_APC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
RCC_APC <- RunHarmony(RCC_APC, group.by.vars = "orig.ident")
RCC_APC <- RunUMAP(RCC_APC, reduction = "harmony", dims = 1:20)
RCC_APC <- RunTSNE(RCC_APC, reduction = "harmony", dims = 1:20)
names(RCC_APC@reductions)

###Ranking.
Idents(RCC_APC) <- RCC_APC$celltype
RCC_APC <- RenameIdents(RCC_APC, "apCAFs" = "apCAFs", "NFs" = "NFs", "Endothelial cells" = "Endothelial cells", "cDC1" = "cDC1", "cDC2" = "cDC2", "C1QC+ Macrophages" = "C1QC+ Macrophages", "SPP1+ Macrophages" = "SPP1+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(RCC_APC))
RCC_APC$celltype <- Idents(RCC_APC)
table(RCC_APC$celltype)
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
saveRDS(RCC_APC, file = "RCC_APC.rds")

###Transcriptomic similarity analysis.
#Returns averaged expression values for each APC cell type.
av <- AverageExpression(RCC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

#Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))

#Using cor() Function to Calculate Spearman Correlation.
cor <- cor(av[cg, ], method = "spearman")

###Visualization.
#Heat map of correlations for Supplementary Figure 14E.
pdf(file="corrplot_RCC_APC_celltype_5000.pdf",width=6, height=6)
corrplot(cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
dev.off()

#2 Supplementary Figure 14F, G, H.
#RCC.
###Loading the Seurat object containing apCAFs and NFs.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
NF_ap <- readRDS(file = "RCC_NF_ap.rds") #The codes for RCC_NF_ap.rds is located at line 172 of Supplementary Figure 9F-G.R.
table(NF_ap$celltype)

###Convert the Seurat object into a CellDataSet object.
Idents(NF_ap) <- NF_ap$celltype
table(Idents(NF_ap))
pd <- NF_ap@meta.data
head(pd)
fd <- data.frame(
  gene_short_name = rownames(NF_ap@assays$RNA) , 
  row.names =  rownames(NF_ap@assays$RNA) 
)
head(fd)
pd <- new("AnnotatedDataFrame",
          data=pd)
fd <- new("AnnotatedDataFrame",
          data=fd)
count=as.data.frame(NF_ap@assays$RNA@counts)
NF_ap_cds <- newCellDataSet(
  as.matrix(count), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=0.5)

###Estimate size factors and dispersions.
NF_ap_cds <- estimateSizeFactors(NF_ap_cds)
NF_ap_cds <- estimateDispersions(NF_ap_cds)

###Filtering low-quality cells.
NF_ap_cds <- detectGenes(NF_ap_cds, min_expr = 0.1) 
print(head(fData(NF_ap_cds)))
NF_ap_cds
expressed_genes <- row.names(subset(fData(NF_ap_cds),
                                    num_cells_expressed >= 10))
print(head(pData(NF_ap_cds)))

###Clustering cells without marker genes.
disp_table <- dispersionTable(NF_ap_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
NF_ap_cds <- setOrderingFilter(NF_ap_cds, unsup_clustering_genes$gene_id)
pdf(file="ordering_genes.pdf",width=5,height=4)
plot_ordering_genes(NF_ap_cds)
dev.off()
pdf(file="pc_variance.pdf",width=8,height=6)
plot_pc_variance_explained(NF_ap_cds, return_all = F) # norm_method='log'
dev.off()
NF_ap_cds <- reduceDimension(NF_ap_cds, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
NF_ap_cds <- clusterCells(NF_ap_cds, num_clusters = 2)
pdf(file="cell_clusters.pdf",width=6,height=5)
plot_cell_clusters(NF_ap_cds, 1, 2, color = "celltype")
dev.off()

###Constructing Single Cell Trajectories.
##1 Choose genes that define a cell's progress.
diff_test_res <- differentialGeneTest(NF_ap_cds[expressed_genes,],
                                      fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
NF_ap_cds <- setOrderingFilter(NF_ap_cds, ordering_genes)
pdf(file="ordering_genes2.pdf",width=5,height=4)
plot_ordering_genes(NF_ap_cds)
dev.off()

##2 Reduce data dimensionality.
NF_ap_cds <- reduceDimension(NF_ap_cds, max_components = 2,
                             method = 'DDRTree')

##3 Order cells along the trajectory.
NF_ap_cds <- orderCells(NF_ap_cds)
saveRDS(NF_ap_cds, file = "NF_ap_cds_RCC.rds")
NF_ap_cds <- orderCells(NF_ap_cds, root_state = 3)

###Visualization.
#plot_cell_trajectory for Supplementary Figure 14F: left.
pdf(file="plot_cell_trajectory_RCC_NF_ap_celltype.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "celltype")  + scale_color_nejm() 
dev.off()

#plot_cell_trajectory for Supplementary Figure 14F: right.
pdf(file="plot_cell_trajectory_RCC_NF_ap_Pseudotime.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "Pseudotime")  
dev.off()

pdf(file="plot_cell_trajectory_RCC_NF_ap_State.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "State")  + scale_color_npg()
dev.off()

###calculate the frequency of distributed cells across different cell subtypes along the pseudotime axis.
mono.info <- pData(NF_ap_cds)
head(mono.info)

#density plot for Supplementary Figure 14G.
pdf(file="plot_density_plot_RCC_NF_ap.pdf",width=5, height=3)
ggplot(mono.info,aes(x=Pseudotime, group=celltype, fill=celltype))+
  geom_density(alpha=0.3,adjust=1.5)+
  scale_fill_manual(values = c(pal_npg("nrc")(10)) )+
  theme_bw()
dev.off()

###Use the function plot_genes_in_pseudotime to draw the expression level of Marker genes.
NF_ap_expressed_genes <-  row.names(subset(fData(NF_ap_cds),
                                           num_cells_expressed >= 10))
NF_ap_filtered <- NF_ap_cds[NF_ap_expressed_genes,]
my_genes <- row.names(subset(fData(NF_ap_filtered),
                             gene_short_name %in% c("CD74")))
cds_subset <- NF_ap_filtered[my_genes,]

#plot_genes_in_pseudotime for Supplementary Figure 14H.
pdf(file="plot_cell_trajectory_RCC_NF_ap_gene.pdf",width=5, height=3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
dev.off()
