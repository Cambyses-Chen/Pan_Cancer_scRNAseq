###codes for Supplementary Figure 9F-G.

library(Seurat)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(dplyr)
library(msigdbr)

#1 Supplementary Figure 9F.
#NPC.
###Load annotated single-cell RNA sequencing data of fibroblasts subpopulations in NPC.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC_Fib_new <- readRDS(file = "NPC_Fib_new.rds") #The codes for NPC_Fib_new.rds is located at line 58 of Supplementary Figure 9A-E.R.
table(NPC_Fib_new$celltype)

###Retain NFs and apCAFs subpopulations.
NF_ap <- subset(NPC_Fib_new, subset = celltype %in% c("NFs", "apCAFs"))
NF_ap$celltype <- as.factor(as.character(NF_ap$celltype))
table(NF_ap$celltype)

###Ranking.
Idents(NF_ap) <- NF_ap$celltype
NF_ap <- RenameIdents(NF_ap, "NFs" = "NFs", "apCAFs" = "apCAFs")
table(Idents(NF_ap))
NF_ap$celltype <- Idents(NF_ap)
table(NF_ap$celltype)
saveRDS(NF_ap, file = "NPC_NF_ap.rds")

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
write.csv(gs.exp_h, file = "gs.exp_NPC_NFs_apCAFs_H.csv")
write.csv(gs.exp_C5GOBP, file = "gs.exp_NPC_NFs_apCAFs_C5GOBP.csv")
saveRDS(gs.exp_h, file = "gs.exp_NPC_NFs_apCAFs_H.rds")
saveRDS(gs.exp_C5GOBP, file = "gs.exp_NPC_NFs_apCAFs_C5GOBP.rds")

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
saveRDS(all, file = "all_H_C5GOBP_NPC_NFs_apCAFs.rds")
write.csv(all, file="all_H_C5GOBP_NPC_NFs_apCAFs.csv",quote=F)

###Choose pathways.
setwd("E:/Raw Data/Choose pathways/NPC")
gs = lapply(readLines("selected_H_GOBP_NPC.txt"), function(x){
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
#Bar plot for Supplementary Figure 9F.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
p <- ggplot(data=dat_plot, aes(x=id, y=t)) + geom_bar(stat="identity") + scale_x_discrete(labels = function(x) str_wrap(x, width = 42))
pdf(file = "NPC_H_GOBP_NFs_apCAFs_barplot.pdf", width=8, height=5)
p + coord_flip() + geom_bar(stat="identity", fill="steelblue") + theme_classic() + theme(axis.title.x = element_text(size = 10)) +
  labs(y = "t value of GSVA score, apCAFs versus NFs", x = "") + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) +
  theme(axis.ticks=element_line(size=1)) + theme(axis.line=element_line(size=1))
dev.off()

#2 Supplementary Figure 9G.
#RCC.
#Load single-cell RNA sequencing data of Fibroblasts major cell type in RCC.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC_Fib <- readRDS(file = "RCC_Fib.rds")  #The codes for RCC_Fib.rds is located at line 316 of Supplementary Figure 5M-N.R.
table(RCC_Fib$Sample)
table(RCC_Fib$Celltype)

#Obtain single-cell RNA sequencing data of normal fibroblasts related to RCC.
RCC_NFs <- subset(RCC_Fib, subset = Sample %in% c("Normal kidney"))
RCC_NFs$Sample <- as.factor(as.character(RCC_NFs$Sample))
RCC_NFs$orig.ident <- as.factor(as.character(RCC_NFs$orig.ident))
table(RCC_NFs$Sample)
table(RCC_NFs$Celltype)
RCC_NFs$celltype <- "NFs"
table(RCC_NFs$celltype)

#Load annotated single-cell RNA sequencing data of CAFs subpopulations in RCC.
RCC_CAFs <- readRDS(file = "RCC_Fib_T_res0.3_annotation_sub.rds")  #The codes for RCC_Fib_T_res0.3_annotation_sub.rds is located at line 436 of Supplementary Figure 5M-N.R.

#Merge into new annotated single-cell RNA sequencing data of fibroblast subpopulations in RCC.
RCC_Fib_new <- merge(x = RCC_CAFs, y = RCC_NFs)

#Ranking.
Idents(RCC_Fib_new) <- RCC_Fib_new$celltype
RCC_Fib_new <- RenameIdents(RCC_Fib_new, "apCAFs" = "apCAFs", "myCAFs" = "myCAFs", "eCAFs" = "eCAFs", "NFs" = "NFs")
table(Idents(RCC_Fib_new))
RCC_Fib_new$celltype <- Idents(RCC_Fib_new)
table(RCC_Fib_new$celltype)

###Normalizing the data.
RCC_Fib_new <- NormalizeData(RCC_Fib_new, normalization.method = "LogNormalize", scale.factor = 10000)

###Identification of highly variable features (feature selection).
RCC_Fib_new <- FindVariableFeatures(RCC_Fib_new, selection.method = "vst", nfeatures = 2000)

###Scaling the data.
all.genes <- rownames(RCC_Fib_new)
RCC_Fib_new <- ScaleData(RCC_Fib_new, features = all.genes)
RCC_Fib_new <- ScaleData(RCC_Fib_new, vars.to.regress = "percent.mt")

###Perform linear dimensional reduction.
RCC_Fib_new <- RunPCA(RCC_Fib_new, features = VariableFeatures(object = RCC_Fib_new))

###harmony integrates data.
###Run non-linear dimensional reduction (UMAP/tSNE).
library(harmony)
RCC_Fib_new <- RunHarmony(RCC_Fib_new, group.by.vars = "orig.ident")
RCC_Fib_new <- RunUMAP(RCC_Fib_new, reduction = "harmony", dims = 1:20)
RCC_Fib_new <- RunTSNE(RCC_Fib_new, reduction = "harmony", dims = 1:20)
names(RCC_Fib_new@reductions)
saveRDS(RCC_Fib_new, file = "RCC_Fib_new.rds")

###Retain NFs and apCAFs subpopulations.
NF_ap <- subset(RCC_Fib_new, subset = celltype %in% c("NFs", "apCAFs"))
NF_ap$celltype <- as.factor(as.character(NF_ap$celltype))
table(NF_ap$celltype)

###Ranking.
Idents(NF_ap) <- NF_ap$celltype
NF_ap <- RenameIdents(NF_ap, "NFs" = "NFs", "apCAFs" = "apCAFs")
table(Idents(NF_ap))
NF_ap$celltype <- Idents(NF_ap)
table(NF_ap$celltype)
saveRDS(NF_ap, file = "RCC_NF_ap.rds")

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
write.csv(gs.exp_h, file = "gs.exp_RCC_NFs_apCAFs_H.csv")
write.csv(gs.exp_C5GOBP, file = "gs.exp_RCC_NFs_apCAFs_C5GOBP.csv")
saveRDS(gs.exp_h, file = "gs.exp_RCC_NFs_apCAFs_H.rds")
saveRDS(gs.exp_C5GOBP, file = "gs.exp_RCC_NFs_apCAFs_C5GOBP.rds")

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
saveRDS(all, file = "all_H_C5GOBP_RCC_NFs_apCAFs.rds")
write.csv(all, file="all_H_C5GOBP_RCC_NFs_apCAFs.csv",quote=F)

###Choose pathways.
setwd("E:/Raw Data/Choose pathways/RCC")
gs = lapply(readLines("selected_H_GOBP_RCC.txt"), function(x){
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
#Bar plot for Supplementary Figure 9G.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
p <- ggplot(data=dat_plot, aes(x=id, y=t)) + geom_bar(stat="identity") + scale_x_discrete(labels = function(x) str_wrap(x, width = 42))
pdf(file = "RCC_H_GOBP_NFs_apCAFs_barplot.pdf", width=8, height=5)
p + coord_flip() + geom_bar(stat="identity", fill="steelblue") + theme_classic() + theme(axis.title.x = element_text(size = 10)) +
  labs(y = "t value of GSVA score, apCAFs versus NFs", x = "") + theme(axis.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) +
  theme(axis.ticks=element_line(size=1)) + theme(axis.line=element_line(size=1))
dev.off()
