###barplot for Figure 1B.
library(ggplot2)
setwd("E:/Raw Data/cancers_information")
t <- read.csv(file = "cancer_sample.CSV")
pdf(file = "cancer_barplot.pdf", width=8, height=6)
ggplot(data=t, aes(x=reorder(Cancer, -Sample), y=Sample)) +
  geom_bar(stat="identity", fill="steelblue") +
  geom_text(aes(label=Sample), vjust=1.6, color="white", size=5)+
  theme_classic() + labs(y = "Number of samples", x = "") + theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1)) + theme(axis.text.y = element_text(size = 20)) +
  theme(axis.ticks=element_line(size=1)) + theme(axis.line=element_line(size=1))
dev.off()
