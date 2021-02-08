library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
#library(Rsamtools)

args <- commandArgs( trailingOnly = TRUE )
reads_file=args[1]
out_file=args[2]

samples=unlist(strsplit(args[3],","))

#!!!TODO: update this and the other qc plot scripts to segregate by batch and by type and include metasheet info

theme_set(theme_classic())
tab <- read.table(reads_file,header=FALSE,sep="\t")
colnames(tab)=c("sample","reads")

tab = subset(tab, tab$sample %in% samples)

ggplot(tab, aes(x =sample,y = reads)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Unique reads mapped to main chromosomes") +
  xlab("Sample")
ggsave(out_file,width = 12, height = 5)
