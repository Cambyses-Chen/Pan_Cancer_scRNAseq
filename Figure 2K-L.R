###codes for Figure 2K-L.

library(Seurat)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(msigdbr)

#1 Figure 2K.
#HNSCC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in HNSCC.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC_Fib_new <- readRDS(file = "HNSCC_Fib_new.rds") #The codes for HNSCC_Fib_new.rds is located at line 58 of Figure 2I-J.R.
table(HNSCC_Fib_new$celltype)

###Retain NFs and apCAFs subpopulations.
NF_ap <- subset(HNSCC_Fib_new, subset = celltype %in% c("NFs", "apCAFs"))
NF_ap$celltype <- as.factor(as.character(NF_ap$celltype))
table(NF_ap$celltype)

###Ranking.
Idents(NF_ap) <- NF_ap$celltype
NF_ap <- RenameIdents(NF_ap, "NFs" = "NFs", "apCAFs" = "apCAFs")
table(Idents(NF_ap))
NF_ap$celltype <- Idents(NF_ap)
table(NF_ap$celltype)
saveRDS(NF_ap, file = "HNSCC_NF_ap.rds")

###Extraction the normalization matrix.
expr <- as.matrix(NF_ap@assays$RNA@data)

###Load gene signature from msigdbr R package.
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
h_gene_sets_list = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
C5GOBP_gene_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)
C5GOBP_gene_sets_list = C5GOBP_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

###Calculate GSVA scores for gene signatures in each cell.
gs.exp_h <- gsva(expr, h_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
gs.exp_C5GOBP <- gsva(expr, C5GOBP_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
write.csv(gs.exp_h, file = "gs.exp_HNSCC_NFs_apCAFs_H.csv")
write.csv(gs.exp_C5GOBP, file = "gs.exp_HNSCC_NFs_apCAFs_C5GOBP.csv")
saveRDS(gs.exp_h, file = "gs.exp_HNSCC_NFs_apCAFs_H.rds")
saveRDS(gs.exp_C5GOBP, file = "gs.exp_HNSCC_NFs_apCAFs_C5GOBP.rds")

###Merge matrices.
identical(colnames(gs.exp_h), colnames(gs.exp_C5GOBP))
gs.exp <- rbind(gs.exp_h, gs.exp_C5GOBP)

###Construct design matrix.
t <- NF_ap@meta.data
t2 <- t[grep("apCAFs", t[,10]),]
t3 <- t[grep("NFs", t[,10]),]
t2$apCAFs <- 1
t2$NFs <- 0
t3$apCAFs <- 0
t3$NFs <- 1
t4 <- rbind(t2, t3)
t5 <- t4[, c(11, 12)]
design <- as.matrix(t5)

###Differential analysis using limma.
library(limma)
compare <- makeContrasts(apCAFs - NFs, levels=design)
fit <- lmFit(gs.exp, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
all <- topTable(fit3, coef=1, number=Inf, adjust.method="BH")
head(all)
saveRDS(all, file = "all_H_C5GOBP_HNSCC_NFs_apCAFs.rds")
write.csv(all, file="all_H_C5GOBP_HNSCC_NFs_apCAFs.csv",quote=F)

###Choose pathways.
setwd("E:/Raw Data/Choose pathways/HNSCC")
gs = lapply(readLines("selected_H_GOBP_HNSCC.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Refine the limma output by selected pathways.
s <- all[gs[[1]],]
dat_plot <- data.frame(id = row.names(s),
                       t = s$t)

###Remove prefixes and underscores from pathway names.
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
dat_plot$id <- str_replace(dat_plot$id , "GOBP_","")
dat_plot$id <- str_replace_all(dat_plot$id, "_", " ")

###Ranking.
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

###Visualization.
#Bar plot for Figure 2K.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
p <- ggplot(data=dat_plot, aes(x=id, y=t)) + geom_bar(stat="identity") + scale_x_discrete(labels = function(x) str_wrap(x, width = 42))
pdf(file = "HNSCC_H_GOBP_NFs_apCAFs_barplot.pdf", width=8, height=5)
p + coord_flip() + geom_bar(stat="identity", fill="steelblue") + theme_classic() + theme(axis.title.x = element_text(size = 10)) +
  labs(y = "t value of GSVA score, apCAFs versus NFs", x = "") + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) +
  theme(axis.ticks=element_line(size=1)) + theme(axis.line=element_line(size=1))
dev.off()

#2 Figure 2L.
#CC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in CC.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CC_Fib_new <- readRDS(file = "CC_Fib_new.rds") #The codes for CC_Fib_new.rds is located at line 118 of Supplementary Figure 9A-E.R.
table(CC_Fib_new$celltype)

###Retain NFs and apCAFs subpopulations.
NF_ap <- subset(CC_Fib_new, subset = celltype %in% c("NFs", "apCAFs"))
NF_ap$celltype <- as.factor(as.character(NF_ap$celltype))
table(NF_ap$celltype)

###Ranking.
Idents(NF_ap) <- NF_ap$celltype
NF_ap <- RenameIdents(NF_ap, "NFs" = "NFs", "apCAFs" = "apCAFs")
table(Idents(NF_ap))
NF_ap$celltype <- Idents(NF_ap)
table(NF_ap$celltype)
saveRDS(NF_ap, file = "CC_NF_ap.rds")

###Extraction the normalization matrix.
expr <- as.matrix(NF_ap@assays$RNA@data)

###Load gene signature from msigdbr R package.
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
h_gene_sets_list = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
C5GOBP_gene_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)
C5GOBP_gene_sets_list = C5GOBP_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

###Calculate GSVA scores for gene signatures in each cell.
gs.exp_h <- gsva(expr, h_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
gs.exp_C5GOBP <- gsva(expr, C5GOBP_gene_sets_list, kcdf="Gaussian",method = "gsva", parallel.sz = 10)
write.csv(gs.exp_h, file = "gs.exp_CC_NFs_apCAFs_H.csv")
write.csv(gs.exp_C5GOBP, file = "gs.exp_CC_NFs_apCAFs_C5GOBP.csv")
saveRDS(gs.exp_h, file = "gs.exp_CC_NFs_apCAFs_H.rds")
saveRDS(gs.exp_C5GOBP, file = "gs.exp_CC_NFs_apCAFs_C5GOBP.rds")

###Merge matrices.
identical(colnames(gs.exp_h), colnames(gs.exp_C5GOBP))
gs.exp <- rbind(gs.exp_h, gs.exp_C5GOBP)

###Construct design matrix.
t <- NF_ap@meta.data
t2 <- t[grep("apCAFs", t[,10]),]
t3 <- t[grep("NFs", t[,10]),]
t2$apCAFs <- 1
t2$NFs <- 0
t3$apCAFs <- 0
t3$NFs <- 1
t4 <- rbind(t2, t3)
t5 <- t4[, c(11, 12)]
design <- as.matrix(t5)

###Differential analysis using limma.
library(limma)
compare <- makeContrasts(apCAFs - NFs, levels=design)
fit <- lmFit(gs.exp, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
all <- topTable(fit3, coef=1, number=Inf, adjust.method="BH")
head(all)
saveRDS(all, file = "all_H_C5GOBP_CC_NFs_apCAFs.rds")
write.csv(all, file="all_H_C5GOBP_CC_NFs_apCAFs.csv",quote=F)

###Choose pathways.
setwd("E:/Raw Data/Choose pathways/CC")
gs = lapply(readLines("selected_H_GOBP_CC.txt"), function(x){
  y = strsplit(x, "\t")[[1]]
  y = y[2:length(y)]
  return(y)
})
gs

###Refine the limma output by selected pathways.
s <- all[gs[[1]],]
dat_plot <- data.frame(id = row.names(s),
                       t = s$t)

###Remove prefixes and underscores from pathway names.
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
dat_plot$id <- str_replace(dat_plot$id , "GOBP_","")
dat_plot$id <- str_replace_all(dat_plot$id, "_", " ")

###Ranking.
dat_plot <- dat_plot %>% arrange(t)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

###Visualization.
#Bar plot for Figure 2L.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
p <- ggplot(data=dat_plot, aes(x=id, y=t)) + geom_bar(stat="identity") + scale_x_discrete(labels = function(x) str_wrap(x, width = 42))
pdf(file = "CC_H_GOBP_NFs_apCAFs_barplot.pdf", width=8, height=5)
p + coord_flip() + geom_bar(stat="identity", fill="steelblue") + theme_classic() + theme(axis.title.x = element_text(size = 10)) +
  labs(y = "t value of GSVA score, apCAFs versus NFs", x = "") + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) +
  theme(axis.ticks=element_line(size=1)) + theme(axis.line=element_line(size=1))
dev.off()
