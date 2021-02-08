args <- commandArgs( trailingOnly = TRUE )
filelist=args[1]
metasheet=args[2]
rms_dir=args[3]

ctdna.file="misc/ctDNA.tsv"

suppressMessages(library(MEDIPS))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(annotatr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(ggdendro))
suppressMessages(library(matrixStats))
suppressMessages(library(stringr))
suppressMessages(library(circlize))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
library(plotROC)

source("scripts/limmamedips.R")

#TODO: review which libraries can be removed

meth.file=file.path(rms_dir, "medip.rds")
rms.file=file.path(rms_dir, "rms.rds") #file will have relative methylation score for all windows in final analysis
restrict=file.path(rms_dir, "reference_up_filtered.bed")
rms.thresh=0.2 #relative methylation score must be greater than this to be counted as methylated

chrs =  paste0("chr", c(1:22, "X", "Y", "M"))
        
#read metasheet	
met=read.table(metasheet, header=T, sep=",", as.is=T)
if(!all(c("SampleName", "Class", "Type") %in% colnames(met))) stop("check metasheet file for required columns")

#will batch info be included?
if("Batch" %in% colnames(met)) show.batch=T

#match files to samples in metasheet
file=unlist(strsplit(filelist," "))
SampleName = file %>% str_replace(".*/","") %>% str_replace(".medip.rds","") 
                                                                                                                          
message("excluding the following samples that are in the config file but not the metasheet file (if any):")       
message(SampleName[!SampleName %in% met$SampleName])  
        
met = merge(met, data.frame(cbind(SampleName, file)))                                                                                  
met$file = as.character(met$file)                                                                                                      
# load medips objects:                                                                                           
message("reading in medips files") 
if(!file.exists(meth.file)) {
obj1 <- lapply(met$file, readRDS)                                                                                

# the coupling set contains counts of CG for each window, can be used for normalization
if(file.exists("analysis/couplingset/couplingset.rds")) {
  CS = readRDS("analysis/couplingset/couplingset.rds")
} else {
  CS = MEDIPS.couplingVector(pattern = "CG", refObj = obj[[1]]) 
}

meth = MEDIPS.meth(obj1, CSet = CS, chr = chrs, MeDIP = T) 
	saveRDS(meth, file = meth.file)
} else {
	meth=readRDS(meth.file)
}

meth.all = meth[,grepl("bam.rms", colnames(meth))]
meth.coords = meth[,c("chr", "start", "stop")]

rms.tot=colSums(meth.all)
rms.tot.binary=colSums(meth.all>rms.thresh)

# restrict to specified sites (ie NEPC-up in LuCaPs):
sites_up=import(restrict, format="BED") 
incl = rep(TRUE, nrow(meth))
incl = incl & (makeGRangesFromDataFrame(meth) %over% sites_up) 	

meth.target = meth[incl,grepl("bam.rms", colnames(meth))]
colnames(meth.target) = str_replace(colnames(meth.target),".dedup.bam.rms","") %>% str_replace("\\.","-")

meth.coords=meth.coords[incl,]

saveRDS(meth.target, rms.file)
score = colSums(meth.target)/rms.tot
score.binary = colSums(meth.target>rms.thresh)/rms.tot.binary

df = data.frame(NE_meth_score = score, NE_meth_score_binary = score.binary, SampleName=names(score))
df = merge(df, met, by="SampleName")
df$Type=as.factor(df$Type)

cols <- c("PRAD" = "blue", "NEPC" = "red")

ggplot(df, aes(y=NE_meth_score, x=Type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=Type), height = 0) + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score.pdf"), height=3, width=3)

if(show.batch) {
  ggplot(df, aes(y=NE_meth_score, x=Type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=as.factor(Batch)), height = 0) + theme_classic()
  ggsave(file.path(rms_dir, "rms_score_batch.pdf"), height=3, width=3)
}

df$State = factor(df$State, levels=c("Localized PRAD", "Metastatic PRAD", "Treatment emergent", "De novo"))
ggplot(df, aes(y=NE_meth_score, x=State)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=State), height = 0) + theme_classic()
ggsave(file.path(rms_dir, "rms_score_state.pdf"), height=3, width=3)

p=wilcox.test(x=df$NE_meth_score[df$Type=="NEPC"], y=df$NE_meth_score[df$Type=="PRAD"])$p.value
message("Mann Whitney test for NEPC vs PRAD NE_meth_scores: ", p)

# plot auc for off-the shelf approach:
curves <-  calculate_roc(D=as.numeric(df$Type=="NEPC"), M=df$NE_meth_score, ci = FALSE)
auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=df, aes(d = as.numeric(df$Type=="NEPC"), m = NE_meth_score), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve.pdf"),
  width = 3, height = 3)

message("classification AUC: ", calc_auc(auc_plot)$AUC)

#write table with results
write.table(df, file = file.path(rms_dir, "rms_scores.txt"), quote=F, sep="\t", row.names=F)

#extra: annotate with ctDNA content estimates (high vs low)
ctdna=read.table(ctdna.file, sep="\t", head=T)
m=merge(df,ctdna, by="SampleName")

#print ctDNA stats
ct = m$ctDNA_estimate
message("cohort ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))
ct = m$ctDNA_estimate[m$Type=="PRAD"]
message("PRAD ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))
ct = m$ctDNA_estimate[m$Type=="NEPC"]
message("NEPC ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))


m$ctDNA = factor(ifelse(m$ctDNA_estimate >= 0.1, ">= 10% ctDNA", " < 10% ctDNA"), levels=c(">= 10% ctDNA", " < 10% ctDNA", "unknown"))
m$ctDNA[is.na(m$ctDNA)] = "unknown"

ggplot(m, aes(y=NE_meth_score, x=Type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=ctDNA), height = 0) + theme_classic()
ggsave(file.path(rms_dir, "rms_score_ctDNA.pdf"), height=3, width=4)

#plot NE score vs ctDNA estimate:
ggplot(m) + geom_point(aes(y=NE_meth_score, x=ctDNA_estimate, col=Type)) + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_ctDNA_scatter.pdf"), height=3, width=4)

#plot NE score vs ctDNA estimate for PRAD only:
ggplot(subset(m, Type=="PRAD"), aes(y=NE_meth_score, x=ctDNA_estimate, col=Type)) + geom_smooth(method="lm", se=F, col="darkgray") + geom_point() + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_ctDNA_scatter_PRAD.pdf"), height=3, width=4)

cor=cor.test(x=m$NE_meth_score[m$Type=="PRAD"], y=m$ctDNA_estimate[m$Type=="PRAD"])
message("correlation for PRAD sample NE score and ctDNA = ", cor$estimate, "(p = ", cor$p.value, ")")

#plot NE score vs ctDNA estimate for NEPC only:
ggplot(subset(m, Type=="NEPC"), aes(y=NE_meth_score, x=ctDNA_estimate, col=Type)) + geom_smooth(method="lm", se=F, col="darkgray") + geom_point(
) + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_ctDNA_scatter_NEPC.pdf"), height=3, width=4)

cor=cor.test(x=m$NE_meth_score[m$Type=="NEPC"], y=m$ctDNA_estimate[m$Type=="NEPC"])
message("correlation for NEPC sample NE score and ctDNA = ", cor$estimate, "(p = ", cor$p.value, ")")
 

#calculate AUC for samples with > 5% ctDNA only
ms = subset(m, ctDNA_estimate > 0.05)

curves <-  calculate_roc(D=as.numeric(ms$Type=="NEPC"), M=ms$NE_meth_score, ci = FALSE)

auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=ms, aes(d = as.numeric(ms$Type=="NEPC"), m = NE_meth_score), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_ctDNA_0.05.pdf"),
  width = 3, height = 3)

ggplot(ms, aes(y=NE_meth_score, x=Type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=Type), height = 0) + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_ctDNA_0.05.pdf"), height=3, width=3)

message("classification AUC for samples with >= 5% ctDNA: ", calc_auc(auc_plot)$AUC)

message("median NE methylation score for PRAD when limiting to >=5% ctDNA: ", median(ms$NE_meth_score[ms$Type=="PRAD"]), ". sd = ", sd(ms$NE_meth_score[ms$Type=="PRAD"]), ", n = ", sum(ms$Type=="PRAD"))

message("median NE methylation score for NEPC when limiting to >=5% ctDNA: ", median(ms$NE_meth_score[ms$Type=="NEPC"]), ". sd = ", sd(ms$NE_meth_score[ms$Type=="NEPC"]), ", n = ", sum(ms$Type=="NEPC"))

message("wilcoxon p val for PRAD vs NEPC = ", wilcox.test(x=ms$NE_meth_score[ms$Type=="NEPC"], y=ms$NE_meth_score[ms$Type=="PRAD"])$p.value)

#do the same for samples with > 10% ctDNA only
ms = subset(m, ctDNA_estimate > 0.10)

curves <- calculate_roc(D=as.numeric(ms$Type=="NEPC"), M=ms$NE_meth_score, ci = FALSE)

auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=ms, aes(d = as.numeric(ms$Type=="NEPC"), m = NE_meth_score), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_ctDNA_0.1.pdf"),
  width = 3, height = 3)

ggplot(ms, aes(y=NE_meth_score, x=Type)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=Type), height = 0) + theme_classic() + scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_ctDNA_0.1.pdf"), height=3, width=3)

message("classification AUC for samples with >= 10% ctDNA: ", calc_auc(auc_plot)$AUC)

message("median NE methylation score for PRAD when limiting to >=10% ctDNA: ", median(ms$NE_meth_score[ms$Type=="PRAD"]), ". sd = ", sd(ms$NE_meth_score[ms$Type=="PRAD"]), ", n = ", sum(ms$Type=="PRAD"))

message("median NE methylation score for NEPC when limiting to >=10% ctDNA: ", median(ms$NE_meth_score[ms$Type=="NEPC"]), ". sd = ", sd(ms$NE_meth_score[ms$Type=="NEPC"]), ", n = ", sum(ms$Type=="NEPC"))

message("wilcoxon p val for PRAD vs NEPC = ", wilcox.test(x=ms$NE_meth_score[ms$Type=="NEPC"], y=ms$NE_meth_score[ms$Type=="PRAD"])$p.value)

