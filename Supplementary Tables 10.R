#1 Supplementary Tables 10 HNSCC slice HNSCC201125T10 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/HNSCC")
HNSCC201125T10 <- readRDS(file = "HNSCC201125T10_HNSCC_semla2.rds") #The detailed code for HNSCC201125T10_HNSCC_semla2.rds is at line 72 of the Supplementary Figure 10A-E.R file.
selected_celltypes <- rownames(HNSCC201125T10)
t <- FetchData(HNSCC201125T10, selected_celltypes)
data <- as.data.frame(t)
write.csv(data, file = "HNSCC_slice_HNSCC201125T10_celltype_proportions.csv")

#2 Supplementary Tables 10 HNSCC slice P210325T3 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/HNSCC")
P210325T3 <- readRDS(file = "P210325T3_HNSCC_semla2.rds") #The detailed code for P210325T3_HNSCC_semla2.rds is at line 72 of the Supplementary Figure 10F-J.R file.
selected_celltypes <- rownames(P210325T3)
t <- FetchData(P210325T3, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "HNSCC_slice_P210325T3_celltype_proportions.csv")

#3 Supplementary Tables 10 HNSCC slice P210325T5 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/HNSCC")
P210325T5 <- readRDS(file = "P210325T5_HNSCC_semla2.rds") #The detailed code for P210325T5_HNSCC_semla2.rds is at line 72 of the Figure 3A-E.R file.
selected_celltypes <- rownames(P210325T5)
t <- FetchData(P210325T5, selected_celltypes)
data <- as.data.frame(t)
write.csv(data, file = "HNSCC_slice_P210325T5_celltype_proportions.csv")

#4 Supplementary Tables 10 OV slice GSM6592133 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/OV")
GSM6592133 <- readRDS(file = "GSM6592133_OV_semla2.rds") #The detailed code for GSM6592133_OV_semla2.rds is at line 72 of the Figure 3F-J.R file.
selected_celltypes <- rownames(GSM6592133)
t <- FetchData(GSM6592133, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "OV_slice_GSM6592133_celltype_proportions.csv")

#5 Supplementary Tables 10 OV slice GSM6592135 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/OV")
GSM6592135 <- readRDS(file = "GSM6592135_OV_semla2.rds") #The detailed code for GSM6592135_OV_semla2.rds is at line 72 of the Supplementary Figure 11A-E.R file.
selected_celltypes <- rownames(GSM6592135)
t <- FetchData(GSM6592135, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "OV_slice_GSM6592135_celltype_proportions.csv")

#6 Supplementary Tables 10 OV slice GSM6592137 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/OV")
GSM6592137 <- readRDS(file = "GSM6592137_OV_semla2.rds") #The detailed code for GSM6592137_OV_semla2.rds is at line 72 of the Supplementary Figure 11F-J.R file.
selected_celltypes <- rownames(GSM6592137)
t <- FetchData(GSM6592137, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "OV_slice_GSM6592137_celltype_proportions.csv")

#7 Supplementary Tables 10 BRCA slice CID4465 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/BRCA")
CID4465 <- readRDS(file = "CID4465_BRCA_semla2.rds") #The detailed code for CID4465_BRCA_semla2.rds is at line 73 of the Supplementary Figure 12A-E.R file.
selected_celltypes <- rownames(CID4465)
t <- FetchData(CID4465, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "BRCA_slice_CID4465_celltype_proportions.csv")

#8 Supplementary Tables 10 BRCA slice CID4535 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/BRCA")
CID4535 <- readRDS(file = "CID4535_BRCA_semla2.rds") #The detailed code for CID4535_BRCA_semla2.rds is at line 73 of the Supplementary Figure 12F-J.R file.
selected_celltypes <- rownames(CID4535)
t <- FetchData(CID4535, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "BRCA_slice_CID4535_celltype_proportions.csv")

#9 Supplementary Tables 10 BRCA slice CID44971 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/BRCA")
CID44971 <- readRDS(file = "CID44971_BRCA_semla2.rds") #The detailed code for CID44971_BRCA_semla2.rds is at line 73 of the Supplementary Figure 12K-O.R file.
selected_celltypes <- rownames(CID44971)
t <- FetchData(CID44971, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "BRCA_slice_CID44971_celltype_proportions.csv")

#10 Supplementary Tables 10 CRC slice GSM7089855 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/CRC")
GSM7089855 <- readRDS(file = "GSM7089855_CRC_semla2.rds") #The detailed code for GSM7089855_CRC_semla2.rds is at line 72 of the Supplementary Figure 13A-E.R file.
selected_celltypes <- rownames(GSM7089855)
t <- FetchData(GSM7089855, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "CRC_slice_GSM7089855_celltype_proportions.csv")

#11 Supplementary Tables 10 CRC slice GSM7089857 sheet: inferring cell type proportions from spatial transcriptome data.
setwd("E:/Raw Data/ST_data/CRC")
GSM7089857 <- readRDS(file = "GSM7089857_CRC_semla2.rds") #The detailed code for GSM7089857_CRC_semla2.rds is at line 72 of the Supplementary Figure 13F-J.R file.
selected_celltypes <- rownames(GSM7089857)
t <- FetchData(GSM7089857, selected_celltypes) 
data <- as.data.frame(t)
write.csv(data, file = "CRC_slice_GSM7089857_celltype_proportions.csv")
