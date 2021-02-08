#plot LFC and model coefficients for DMRs
library(ggplot2)
library(dplyr)

args <- commandArgs( trailingOnly = TRUE )

bed=args[1]
out.dir=args[2]

t=read.table(bed, sep = "\t", header=T, stringsAsFactors=F)
t$sig.plasma = as.factor(ifelse(t$sig.plasma==1,"significant", "not signficant"))
t$dmr = paste0(t$chr,":", t$start, "-", t$end)

ggplot(t, aes(y=lfc.plasma, x=lfc.lucap, color=sig.plasma)) + geom_point(pch=".") + theme_classic()
ggsave(file.path(out.dir, "cf_vs_lucap.pdf"), height=4, width=4)

corr = cor.test (~ abs(lfc.plasma) + abs(lfc.lucap), data = t, method="spearman")
message("Spearman correlation coefficient between LuCaP and cfDNA LFC is ", corr$estimate, ". p=", corr$p.value)

#how many times each DMR is included in a model
n = t %>% group_by(dmr) %>% summarize(n=n()) %>% data.frame()
n = n[order(-n$n),]
n$dmr_order = seq(1,nrow(n))

#how many times a DMR is significant
s = t %>% group_by(dmr) %>% summarize(nsig = sum(sig.plasma=="significant")) %>% data.frame()
n = merge(s,n,by="dmr")

ggplot(n) + geom_col(aes(x=dmr_order,y=n)) + geom_point(aes(x=dmr_order,y=nsig), shape="-") + theme_classic()
ggsave(file.path(out.dir,"dmr.hist.pdf"), height=4, width=4)

