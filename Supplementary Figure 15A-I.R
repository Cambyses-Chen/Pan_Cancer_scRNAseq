###codes for Supplementary Figure 15A-I.

#1 Supplementary Figure 15A.
#TCGA-HNSC
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 97 of Figure 2A-D.R.
setwd("E:/Raw Data/cibersort_results/TCGA-HNSC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-HNSC_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/HNSCC/GSE181919")
#scatter plot for Supplementary Figure 15A.
pdf(file = "HNSCC_TCGA-HNSC_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#2 Supplementary Figure 15B.
#TCGA-OV
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 303 of Figure 2A-D.R.
setwd("E:/Raw Data/cibersort_results/TCGA-OV")
ciber <- read.csv(file = "CIBERSORTx_TCGA-OV_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/OV/GSE165897")
#scatter plot for Supplementary Figure 15B.
pdf(file = "OV_TCGA-OV_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#3 Supplementary Figure 15C.
#GSE102349_NPC
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 62 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/GSE102349_NPC_RNAseq")
ciber <- read.csv(file = "CIBERSORTx_GSE102349_NPC_RNAseq_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/NPC/HRA000087")
#scatter plot for Supplementary Figure 15C.
pdf(file = "NPC_GSE102349_NPC_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#4 Supplementary Figure 15D.
#TCGA-CESC
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 232 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-CESC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-CESC_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/CC/GSE208653_RAW")
#scatter plot for Supplementary Figure 15D.
pdf(file = "CC_TCGA-CESC_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#5 Supplementary Figure 15E.
#TCGA-COAD
### Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 443 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-COAD")
ciber <- read.csv(file = "CIBERSORTx_TCGA-COAD_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
#scatter plot for Supplementary Figure 15E.
pdf(file = "CRC_TCGA-COAD_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#6 Supplementary Figure 15F.
#TCGA-READ
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 616 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-READ")
ciber <- read.csv(file = "CIBERSORTx_TCGA-READ_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/CRC/HRA000979")
#scatter plot for Supplementary Figure 15F.
pdf(file = "CRC_TCGA-READ_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#7 Supplementary Figure 15G.
#TCGA-BRCA
###Correlation analysis.
###Load and merge two parts of the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 832 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-BRCA")
ciber_p1 <- read.csv(file = "CIBERSORTx_TCGA-BRCA-p1_Results.csv")
ciber_p2 <- read.csv(file = "CIBERSORTx_TCGA-BRCA-p2_Results.csv")
ciber <- rbind(ciber_p1, ciber_p2)
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/BRCA/GSE176078")
#scatter plot for Supplementary Figure 15G.
pdf(file = "BRCA_TCGA-BRCA_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#8 Supplementary Figure 15H.
#TCGA-SKCM
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 1045 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-SKCM")
ciber <- read.csv(file = "CIBERSORTx_TCGA-SKCM_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/CM/GSE215120")
#scatter plot for Supplementary Figure 15H.
pdf(file = "CM_TCGA-SKCM_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()

#9 Supplementary Figure 15I.
#TCGA-KIRC
###Correlation analysis.
###loading the CIBERSORTx result. The detailed code for CIBERSORTx result is located at line 1255 of Supplementary Figure 7A-N.R.
setwd("E:/Raw Data/cibersort_results/TCGA-KIRC")
ciber <- read.csv(file = "CIBERSORTx_TCGA-KIRC_Results.csv")
library(ggplot2)
library(ggpubr)
data <- ciber
data <- as.data.frame(data)

###Shapiro-Wilk test for normality.
shapiro.test(data$apCAFs)
shapiro.test(data$CD8..T.cells)

###Visualization.
setwd("E:/Raw Data/scRNA-seq data/RCC/EGAD00001008030")
#scatter plot for Supplementary Figure 15I.
pdf(file = "RCC_TCGA-KIRC_cibersortx_apCAFs_CD8T_spearman.pdf", width=6, height=5)
ggscatter(data, x = "CD8..T.cells", y = "apCAFs", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "CD8+ T cells signature score", ylab = "apCAFs signature score", cor.coef.size = 8, color = "blue", size = 3) +
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("xy.text", size = 20)
dev.off()
