library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs( trailingOnly = TRUE )
mqc_file=args[1]
out_file=args[2]
samples=unlist(strsplit(args[3],","))

theme_set(theme_classic())

bowtie <- read_tsv(mqc_file)
bowtie$Sample = gsub(".bowtie2", "", bowtie$Sample)

bowtie=bowtie[bowtie$Sample %in% samples,]
ggplot(bowtie, aes(x = Sample,y = overall_alignment_rate)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Bowtie2 Overall Alignment Rate") +
  xlab("Sample")
ggsave(out_file,width = 15, height = 5)
