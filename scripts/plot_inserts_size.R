library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

cat("plotting insert sizes\n")
args <- commandArgs( trailingOnly = TRUE )
insert_files=unlist(strsplit(args[1]," "))
out_file=args[2]
samples=unlist(strsplit(args[3],","))

#function to get sample names from insert size file names
get_sample_name = function(file_name) {
	return(file_name %>% str_replace(".*/","") %>% str_replace(".insert_sizes.txt",""))
}

theme_set(theme_classic())
df = data.frame(sample=character(0), insert_size=numeric(0))
for (i in 1:length(insert_files)) {
	sn = get_sample_name(insert_files[i])
	if (sn %in% samples) {
		ins = read.table(insert_files[i],sep="\t",header=FALSE)
		df = rbind(df, data.frame(sample=rep(sn, nrow(ins)), insert_size = ins$V1))
	}
}

ggplot(df, aes(x =sample,y = insert_size)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Insert size") +
  xlab("Sample") + ylim(0,1000)
ggsave(out_file,width = 15, height = 5)
