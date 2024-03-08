###codes for Supplementary Figure 14A-D.

library(Seurat)
library(scales)
library(ggplot2)
library(corrplot)
library(ggsci)
library(reshape2)
library(monocle)
library(ggridges)

#1 Supplementary Figure 14A.
#NPC.
###Load Seurat objects of various annotated cell types from NPC single-cell RNA sequencing data.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC <- readRDS(file = "NPC_resolution0.2_annotation_sub.rds") #The codes for NPC_resolution0.2_annotation_sub.rds is located at line 616 of Supplementary Figure 5A-B.R.
NPC_M <- readRDS(file = "NPC_Myeloid_res0.3_annotation_sub.rds") #The codes for NPC_Myeloid_res0.3_annotation_sub.rds is located at line 187 of Supplementary Figure 6E-F.R.
NPC_NF_ap <- readRDS(file = "NPC_NF_ap.rds") #The codes for NPC_NF_ap.rds is located at line 28 of Supplementary Figure 9F-G.R.
table(NPC$Celltype)
table(NPC_M$celltype)
table(NPC_NF_ap$celltype)

###Extract endothelial cells.
NPC_Endo <- subset(NPC, subset=Celltype %in% c("Endothelial cells"))
NPC_Endo$Celltype <- as.factor(as.character(NPC_Endo$Celltype))
table(NPC_Endo$Celltype)
NPC_Endo$celltype <- NPC_Endo$Celltype

###Extract the required myeloid cells.
NPC_M2 <- subset(NPC_M, subset=celltype %in% c("cDC1", "cDC2", "LAMP3+ DCs", "C1QC+ Macrophages", "VCAN+ Macrophages"))
NPC_M2$celltype <- as.factor(as.character(NPC_M2$celltype))
table(NPC_M2$celltype)

###Merge into a Seurat object containing various antigen-presenting cells (APCs).
NPC_APC <- merge(x = NPC_Endo, y = c(NPC_M2, NPC_NF_ap))
table(NPC_APC$celltype)

###Normalizing the data.
NPC_APC <- NormalizeData(NPC_APC, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
NPC_APC <- FindVariableFeatures(NPC_APC, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(NPC_APC)
NPC_APC <- ScaleData(NPC_APC, features = all.genes)
NPC_APC <- ScaleData(NPC_APC, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
NPC_APC <- RunPCA(NPC_APC, features = VariableFeatures(object = NPC_APC))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
NPC_APC <- RunHarmony(NPC_APC, group.by.vars = "orig.ident")
NPC_APC <- RunUMAP(NPC_APC, reduction = "harmony", dims = 1:20)
NPC_APC <- RunTSNE(NPC_APC, reduction = "harmony", dims = 1:20)
names(NPC_APC@reductions)

###Ranking.
Idents(NPC_APC) <- NPC_APC$celltype
NPC_APC <- RenameIdents(NPC_APC, "apCAFs" = "apCAFs", "NFs" = "NFs", "Endothelial cells" = "Endothelial cells", "cDC1" = "cDC1", "cDC2" = "cDC2", "LAMP3+ DCs" = "LAMP3+ DCs", "C1QC+ Macrophages" = "C1QC+ Macrophages", "VCAN+ Macrophages" = "VCAN+ Macrophages")
table(Idents(NPC_APC))
NPC_APC$celltype <- Idents(NPC_APC)
table(NPC_APC$celltype)
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
saveRDS(NPC_APC, file = "NPC_APC.rds")

###Transcriptomic similarity analysis.
#Returns averaged expression values for each APC cell type.
av <- AverageExpression(NPC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

#Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))

#Using cor() Function to Calculate Spearman Correlation.
cor <- cor(av[cg, ], method = "spearman")

###Visualization.
#Heat map of correlations for Supplementary Figure 14A.
pdf(file="corrplot_NPC_APC_celltype_5000.pdf",width=6, height=6)
corrplot(cor, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)
dev.off()

#2 Supplementary Figure 14B, C, D.
#NPC.
###Loading the Seurat object containing apCAFs and NFs.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NF_ap <- readRDS(file = "NPC_NF_ap.rds") #The codes for NPC_NF_ap.rds is located at line 28 of Supplementary Figure 9F-G.R.
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
NF_ap_cds <- reduceDimension(NF_ap_cds, max_components = 2, num_dim = 28,
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
saveRDS(NF_ap_cds, file = "NF_ap_cds_NPC.rds")
NF_ap_cds <- orderCells(NF_ap_cds, root_state = 5)

###Visualization.
#plot_cell_trajectory for Supplementary Figure 14B: left.
pdf(file="plot_cell_trajectory_NPC_NF_ap_celltype.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "celltype")  + scale_color_nejm() 
dev.off()

#plot_cell_trajectory for Supplementary Figure 14B: right.
pdf(file="plot_cell_trajectory_NPC_NF_ap_Pseudotime.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "Pseudotime")  
dev.off()

pdf(file="plot_cell_trajectory_NPC_NF_ap_State.pdf",width=4, height=3)
plot_cell_trajectory(NF_ap_cds, color_by = "State")  + scale_color_npg()
dev.off()

###calculate the frequency of distributed cells across different cell subtypes along the pseudotime axis.
mono.info <- pData(NF_ap_cds)
head(mono.info)

#density plot for Supplementary Figure 14C.
pdf(file="plot_density_plot_NPC_NF_ap.pdf",width=5, height=3)
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

#plot_genes_in_pseudotime for Supplementary Figure 14D.
pdf(file="plot_cell_trajectory_NPC_NF_ap_gene.pdf",width=5, height=3)
plot_genes_in_pseudotime(cds_subset, color_by = "celltype") + scale_color_nejm()
dev.off()
