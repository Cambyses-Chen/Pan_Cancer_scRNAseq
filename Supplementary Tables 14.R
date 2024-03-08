#1 Supplementary Tables 14 HNSCC Fib HALLMARK_GLYCOLYSIS sheet: GSVA results for HALLMARK_GLYCOLYSIS pathway in fibroblasts subpopulations derived from scRNA-seq data.
###GSVA scores for gene signatures in each cell.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
gs.exp_h <- readRDS(file = "gs.exp_HNSCC_Fib_H.rds") #The detailed code for gs.exp_HNSCC_Fib_H.rds is at line 151 of the Figure 6A-F.R file.

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
write.csv(exp_g, file = "exp_g_HALLMARK_GLYCOLYSIS_HNSCC_GSVA.csv")

#2 Supplementary Tables 14 NPC Fib HALLMARK_GLYCOLYSIS sheet: GSVA results for HALLMARK_GLYCOLYSIS pathway in fibroblasts subpopulations derived from scRNA-seq data.
###GSVA scores for gene signatures in each cell.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
gs.exp_h <- readRDS(file = "gs.exp_NPC_Fib_H.rds") #The detailed code for gs.exp_NPC_Fib_H.rds is at line 200 of the Figure 6A-F.R file.

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
write.csv(exp_g, file = "exp_g_HALLMARK_GLYCOLYSIS_NPC_GSVA.csv")

#3 Supplementary Tables 14 BRCA Fib HALLMARK_GLYCOLYSIS sheet: GSVA results for HALLMARK_GLYCOLYSIS pathway in fibroblasts subpopulations derived from scRNA-seq data.
###GSVA scores for gene signatures in each cell.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
gs.exp_h <- readRDS(file = "gs.exp_BRCA_Fib_H.rds") #The detailed code for gs.exp_BRCA_Fib_H.rds is at line 249 of the Figure 6A-F.R file.

###select HALLMARK_GLYCOLYSIS GSVA scores.
exp_g <- t(gs.exp_h["HALLMARK_GLYCOLYSIS", , drop = FALSE])
exp_g <- as.data.frame(exp_g)
write.csv(exp_g, file = "exp_g_HALLMARK_GLYCOLYSIS_BRCA_GSVA.csv")
