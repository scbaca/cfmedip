library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs( trailingOnly = TRUE )
mqc_file=args[1]
out_file=args[2]
samples=unlist(strsplit(args[3],","))

theme_set(theme_classic())

tab <- read_tsv(mqc_file)

bowtie=tab[grepl("bowtie",tab$Sample),]%>%
  mutate(id = gsub(".bowtie2", "", Sample)) %>%
  dplyr::rename(rate = `Bowtie 2_mqc-generalstats-bowtie_2-overall_alignment_rate`)

bowtie=bowtie[bowtie$id %in% samples,]
ggplot(bowtie, aes(x =id,y = rate)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Bowtie2 Overall Alignment Rate") +
  xlab("Sample")
ggsave(out_file,width = 12, height = 5)
