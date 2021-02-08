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
fastq=tab[!grepl("bowtie|_2",tab$Sample),]%>%
  mutate(id = gsub("_1", "", Sample)) %>%
  mutate(Duplicate = `FastQC_mqc-generalstats-fastqc-total_sequences` * 
    `FastQC_mqc-generalstats-fastqc-percent_duplicates`/100,
      Unique = `FastQC_mqc-generalstats-fastqc-total_sequences` - 
       Duplicate) %>% gather(type, count, 10:11)
fastq=fastq[fastq$id %in% samples,]

ggplot(fastq, aes(x = id, y = count, fill=type)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Number of Reads")
ggsave(out_file,width = 12, height = 5)
