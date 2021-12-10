library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs( trailingOnly = TRUE )
CpG_file=args[1]
out_file=args[2]
samples = args[3]

theme_set(theme_classic())
samples = unlist(strsplit(args[3],","))
tab <- read.table(CpG_file,colClasses=c("character","NULL","numeric"),header=FALSE,sep=" ")
colnames(tab)=c("sample","enrichment")

tab = subset(tab, sample %in% samples)

ggplot(tab, aes(x =sample,y = enrichment)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CpG enrichment in reads") +
  xlab("Sample")
ggsave(out_file,width = 15, height = 5)
