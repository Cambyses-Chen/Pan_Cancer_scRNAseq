#1 Supplementary Tables 12 HNSCC Top5000 genes sheet: The top 5000 variably expressed genes in apCAFs, as well as their potential source cell types such as NFs, endothelial cells, different macrophages, and various dendritic cells, were identified using scRNA-seq data.

setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_APC <- readRDS(file = "HNSCC_APC.rds") #The detailed code for HNSCC_APC.rds is at line 62 of the Figure 4A-D.R file.

###Transcriptomic similarity analysis.
###Returns averaged expression values for each APC cell type.
av <- AverageExpression(HNSCC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

###Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))
top5000 <- av[cg, ]
write.csv(top5000, file = "top5000_HNSCC_APC.csv")

#2 Supplementary Tables 12 CC Top5000 genes sheet: The top 5000 variably expressed genes in apCAFs, as well as their potential source cell types such as NFs, endothelial cells, different macrophages, and various dendritic cells, were identified using scRNA-seq data.

setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CC_APC <- readRDS(file = "CC_APC.rds") #The detailed code for CC_APC.rds is at line 67 of the Figure 4E-H.R file.

###Transcriptomic similarity analysis.
###Returns averaged expression values for each APC cell type.
av <- AverageExpression(CC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

###Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))
top5000 <- av[cg, ]
write.csv(top5000, file = "top5000_CC_APC.csv")

#3 Supplementary Tables 12 NPC Top5000 genes sheet: The top 5000 variably expressed genes in apCAFs, as well as their potential source cell types such as NFs, endothelial cells, different macrophages, and various dendritic cells, were identified using scRNA-seq data.

setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_APC <- readRDS(file = "NPC_APC.rds") #The detailed code for NPC_APC.rds is at line 67 of the Supplementary Figure 14A-D.R file.

###Transcriptomic similarity analysis.
###Returns averaged expression values for each APC cell type.
av <- AverageExpression(NPC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

###Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))
top5000 <- av[cg, ]
write.csv(top5000, file = "top5000_NPC_APC.csv")

#4 Supplementary Tables 12 RCC Top5000 genes sheet: The top 5000 variably expressed genes in apCAFs, as well as their potential source cell types such as NFs, endothelial cells, different macrophages, and various dendritic cells, were identified using scRNA-seq data.

setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC_APC <- readRDS(file = "RCC_APC.rds") #The detailed code for RCC_APC.rds is at line 67 of the Supplementary Figure 14E-H.R file.

###Transcriptomic similarity analysis.
###Returns averaged expression values for each APC cell type.
av <- AverageExpression(RCC_APC, group.by = "celltype", assays = "RNA")
av = av[[1]]

###Selecting the top 5000 genes with the highest level of dispersion.
cg = names(tail(sort(apply(av, 1, sd)), 5000))
top5000 <- av[cg, ]
write.csv(top5000, file = "top5000_RCC_APC.csv")
