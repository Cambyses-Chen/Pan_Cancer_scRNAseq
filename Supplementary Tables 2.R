#1 Supplementary Tables 2 HNSCC sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
HNSCC <- readRDS(file = "HNSCC_resolution0.3_annotation_sub.rds") #The detailed code for HNSCC_resolution0.3_annotation_sub.rds is at line 273 of the Figure 1C-F.R file.
table(HNSCC$Celltype)

###DEG in Celltype.
Idents(HNSCC) <- HNSCC$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
HNSCC.markers <- FindAllMarkers(object = HNSCC,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=HNSCC.markers[(abs(as.numeric(as.vector(HNSCC.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(HNSCC.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 HNSCC sheet.
write.csv(sig.markers, file="HNSCC.markers_Celltype.csv")

#2 Supplementary Tables 2 OV sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
OV <- readRDS(file = "OV_resolution0.2_annotation2.rds") #The detailed code for OV_resolution0.2_annotation2.rds is at line 261 of the Figure 1G-J.R file.
table(OV$Celltype)

###DEG in Celltype.
Idents(OV) <- OV$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
OV.markers <- FindAllMarkers(object = OV,
                             only.pos = FALSE,
                             min.pct = 0.25,
                             logfc.threshold = logFCfilter)
sig.markers=OV.markers[(abs(as.numeric(as.vector(OV.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(OV.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 OV sheet.
write.csv(sig.markers, file="OV.markers_Celltype.csv")

#3 Supplementary Tables 2 NPC sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
NPC <- readRDS(file = "NPC_resolution0.2_annotation_sub.rds") #The detailed code for NPC_resolution0.2_annotation_sub.rds is at line 616 of the Supplementary Figure 5A-B.R file.
table(NPC$Celltype)

###DEG in Celltype.
Idents(NPC) <- NPC$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
NPC.markers <- FindAllMarkers(object = NPC,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=NPC.markers[(abs(as.numeric(as.vector(NPC.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(NPC.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 NPC sheet.
write.csv(sig.markers, file="NPC.markers_Celltype.csv")

#4 Supplementary Tables 2 CC sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
CC <- readRDS(file = "Cervical_resolution0.3_annotation_2.rds") #The detailed code for Cervical_resolution0.3_annotation_2.rds is at line 268 of the Supplementary Figure 5C-D.R file.
table(CC$Celltype)

###DEG in Celltype.
Idents(CC) <- CC$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
CC.markers <- FindAllMarkers(object = CC,
                             only.pos = FALSE,
                             min.pct = 0.25,
                             logfc.threshold = logFCfilter)
sig.markers=CC.markers[(abs(as.numeric(as.vector(CC.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CC.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 CC sheet.
write.csv(sig.markers, file="CC.markers_Celltype.csv")

#5 Supplementary Tables 2 CRC sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
CRC <- readRDS(file = "CRC_HRA000979_resolution0.1_annotation_sub.rds") #The detailed code for CRC_HRA000979_resolution0.1_annotation_sub.rds is at line 269 of the Supplementary Figure 5E-F.R file.
table(CRC$Celltype)

###DEG in Celltype.
Idents(CRC) <- CRC$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
CRC.markers <- FindAllMarkers(object = CRC,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=CRC.markers[(abs(as.numeric(as.vector(CRC.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CRC.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 CRC sheet.
write.csv(sig.markers, file="CRC.markers_Celltype.csv")

#6 Supplementary Tables 2 BRCA sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
BRCA <- readRDS(file = "BC_GSE176078_resolution0.1_annotation_sub.rds") #The detailed code for BC_GSE176078_resolution0.1_annotation_sub.rds is at line 308 of the Supplementary Figure 5G-H.R file.
table(BRCA$Celltype)

###DEG in Celltype.
Idents(BRCA) <- BRCA$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
BRCA.markers <- FindAllMarkers(object = BRCA,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=BRCA.markers[(abs(as.numeric(as.vector(BRCA.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(BRCA.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 BRCA sheet.
write.csv(sig.markers, file="BRCA.markers_Celltype.csv")

#7 Supplementary Tables 2 AM sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/AM/GSE189889")
AM <- readRDS(file = "AM_resolution0.1_annotation_sub.rds") #The detailed code for AM_resolution0.1_annotation_sub.rds is at line 436 of the Supplementary Figure 5I-J.R file.
table(AM$Celltype)

###DEG in Celltype.
Idents(AM) <- AM$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
AM.markers <- FindAllMarkers(object = AM,
                             only.pos = FALSE,
                             min.pct = 0.25,
                             logfc.threshold = logFCfilter)
sig.markers=AM.markers[(abs(as.numeric(as.vector(AM.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(AM.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 AM sheet.
write.csv(sig.markers, file="AM.markers_Celltype.csv")

#8 Supplementary Tables 2 CM sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
CM <- readRDS(file = "CM_resolution0.2_annotation2.rds") #The detailed code for CM_resolution0.2_annotation2.rds is at line 246 of the Supplementary Figure 5K-L.R file.
table(CM$Celltype)

###DEG in Celltype.
Idents(CM) <- CM$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
CM.markers <- FindAllMarkers(object = CM,
                             only.pos = FALSE,
                             min.pct = 0.25,
                             logfc.threshold = logFCfilter)
sig.markers=CM.markers[(abs(as.numeric(as.vector(CM.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CM.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 CM sheet.
write.csv(sig.markers, file="CM.markers_Celltype.csv")

#9 Supplementary Tables 2 RCC sheet.
library(Seurat)
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
RCC <- readRDS(file = "RCC_resolution0.1_annotation2.rds") #The detailed code for RCC_resolution0.1_annotation2.rds is at line 288 of the Supplementary Figure 5M-N.R file.
table(RCC$Celltype)

###DEG in Celltype.
Idents(RCC) <- RCC$Celltype
logFCfilter=0.5
adjPvalFilter=0.05
RCC.markers <- FindAllMarkers(object = RCC,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter)
sig.markers=RCC.markers[(abs(as.numeric(as.vector(RCC.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(RCC.markers$p_val_adj))<adjPvalFilter),]

###Save the table information for Supplementary Tables 2 RCC sheet.
write.csv(sig.markers, file="RCC.markers_Celltype.csv")
