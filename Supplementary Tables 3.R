#1 Supplementary Tables 3 HNSCC CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
CAF <- readRDS(file = "HNSCC_Fib_T_res0.2_annotation_sub.rds") #The detailed code for HNSCC_Fib_T_res0.2_annotation_sub.rds is at line 418 of the Figure 1C-F.R file.
table(CAF$Sample)
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 HNSCC CAFs sheet.
write.csv(sig.markers, file="HNSCC_CAFs.markers_celltype.csv")

#2 Supplementary Tables 3 OV CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
CAF <- readRDS(file = "OV_Fib_res0.2_annotation_sub.rds") #The detailed code for OV_Fib_res0.2_annotation_sub.rds is at line 406 of the Figure 1G-J.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 OV CAFs sheet.
write.csv(sig.markers, file="OV_CAFs.markers_celltype.csv")

#3 Supplementary Tables 3 NPC CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
CAF <- readRDS(file = "NPC_Fib_T_res0.1_annotation_sub.rds") #The detailed code for NPC_Fib_T_res0.1_annotation_sub.rds is at line 771 of the Supplementary Figure 5A-B.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 NPC CAFs sheet.
write.csv(sig.markers, file="NPC_CAFs.markers_celltype.csv")

#4 Supplementary Tables 3 CC CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CAF <- readRDS(file = "Cervical_Fib_T_res1_annotation_sub.rds") #The detailed code for Cervical_Fib_T_res1_annotation_sub.rds is at line 418 of the Supplementary Figure 5C-D.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 CC CAFs sheet.
write.csv(sig.markers, file="CC_CAFs.markers_celltype.csv")

#5 Supplementary Tables 3 CRC CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
CAF <- readRDS(file = "CRC_HRA000979_Fib_T_res0.2_annotation.rds") #The detailed code for CRC_HRA000979_Fib_T_res0.2_annotation.rds is at line 411 of the Supplementary Figure 5E-F.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 CRC CAFs sheet.
write.csv(sig.markers, file="CRC_CAFs.markers_celltype.csv")

#6 Supplementary Tables 3 BRCA CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
CAF <- readRDS(file = "BC_GSE176078_Fib_res0.2_annotation2.rds") #The detailed code for BC_GSE176078_Fib_res0.2_annotation2.rds is at line 437 of the Supplementary Figure 5G-H.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 BRCA CAFs sheet.
write.csv(sig.markers, file="BRCA_CAFs.markers_celltype.csv")

#7 Supplementary Tables 3 AM CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
CAF <- readRDS(file = "AM_Fib_res0.1_annotation_sub.rds") #The detailed code for AM_Fib_res0.1_annotation_sub.rds is at line 577 of the Supplementary Figure 5I-J.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 AM CAFs sheet.
write.csv(sig.markers, file="AM_CAFs.markers_celltype.csv")

#8 Supplementary Tables 3 CM CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CAF <- readRDS(file = "CM_Fib_res0.4_annotation_sub.rds") #The detailed code for CM_Fib_res0.4_annotation_sub.rds is at line 392 of the Supplementary Figure 5K-L.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 CM CAFs sheet.
write.csv(sig.markers, file="CM_CAFs.markers_celltype.csv")

#9 Supplementary Tables 3 RCC CAFs sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
CAF <- readRDS(file = "RCC_Fib_T_res0.3_annotation_sub.rds")  #The detailed code for RCC_Fib_T_res0.3_annotation_sub.rds is at line 436 of the Supplementary Figure 5M-N.R file.
table(CAF$celltype)

###DEG in CAFs.
Idents(CAF) <- CAF$celltype
logFCfilter=0.5
adjPvalFilter=0.05
CAF.markers <- FindAllMarkers(object = CAF,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CAF.markers[(abs(as.numeric(as.vector(CAF.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CAF.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 3 RCC CAFs sheet.
write.csv(sig.markers, file="RCC_CAFs.markers_celltype.csv")
