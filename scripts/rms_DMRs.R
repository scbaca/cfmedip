args <- commandArgs( trailingOnly = TRUE )
filelist=args[1]
metasheet=args[2]
rms_dir=args[3]

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
rms.file.up=file.path(rms_dir, "rms.up.rds") #file will have relative methylation score for windows that are UP in reference cases
restrict.up=file.path(rms_dir, "reference_up_filtered.bed")
rms.file.down=file.path(rms_dir, "rms.down.rds") #file will have relative methylation score for windows that are DOWN in reference cases
restrict.down=file.path(rms_dir, "reference_down_filtered.bed")
#rms.thresh=0.2 #relative methylation score must be greater than this to be counted as methylated

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

rms.tot=colSums(meth.all)

# restrict to sites that are UP in reference cases
sites_up=import(restrict.up, format="BED") 
incl.up = rep(TRUE, nrow(meth))
incl.up = incl.up & (makeGRangesFromDataFrame(meth) %over% sites_up) 	

meth.target.up = meth[incl.up,grepl("bam.rms", colnames(meth))]
colnames(meth.target.up) = str_replace(colnames(meth.target.up),".dedup.bam.rms","") %>% str_replace("\\.","-")
saveRDS(meth.target.up, rms.file.up)

# restrict to sites that are DOWN in reference cases
sites_down=import(restrict.down, format="BED")
incl.down = rep(TRUE, nrow(meth))
incl.down = incl.down & (makeGRangesFromDataFrame(meth) %over% sites_down)

meth.target.down = meth[incl.down,grepl("bam.rms", colnames(meth))]
colnames(meth.target.down) = str_replace(colnames(meth.target.down),".dedup.bam.rms","") %>% str_replace("\\.","-")
saveRDS(meth.target.down, rms.file.down)

score.up = colSums(meth.target.up)/rms.tot
score.down = colSums(meth.target.down)/rms.tot

df = data.frame(NE_meth_score = score.up, PRAD_meth_score = score.down, SampleName=names(score.up))
df = merge(df, met, by="SampleName")
df$Type=factor(df$Type, levels=c("HEALTHY","PRAD", "NEPC"))

df$portion=log((df$NE_meth_score/sum(incl.up))/(df$PRAD_meth_score/sum(incl.down)),2)

#normalize portion score to the median for the healthy controls:
med_portion = median(df$portion[df$Type=="HEALTHY"])
df$portion = df$portion - med_portion


#write results to text file
write.table(df, file = file.path(rms_dir, "rms_scores.txt"), quote=F, sep="\t", row.names=F)

#subset to remove samples with <0.03 ctDNA

df = subset(df, ctDNA > 0.03 & Type != "HEALTHY")


cols <- c("PRAD" = "blue", "NEPC" = "red")

ggplot(subset(df, Type != "HEALTHY"), aes(y=NE_meth_score, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_NEPC_all.pdf"), height=3, width=3)

ggplot(subset(df, Type != "HEALTHY"), aes(y=PRAD_meth_score, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_PRAD_all.pdf"), height=3, width=3)

ggplot(subset(df, Type != "HEALTHY"), aes(y=portion, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_ratio_all.pdf"), height=3, width=3)

# -----

message("combined cohorts")
p=wilcox.test(x=df$NE_meth_score[df$Type=="NEPC"], y=df$NE_meth_score[df$Type=="PRAD"])$p.value    
message("Mann Whitney test for NEPC vs PRAD NE_meth_scores: ", p) 

p=wilcox.test(x=df$PRAD_meth_score[df$Type=="NEPC"], y=df$PRAD_meth_score[df$Type=="PRAD"])$p.value                
message("Mann Whitney test for NEPC vs PRAD PRAD_meth_scores: ", p) 

p=wilcox.test(x=df$portion[df$Type=="NEPC"], y=df$portion[df$Type=="PRAD"])$p.value                
message("Mann Whitney test for NEPC vs PRAD portion scores: ", p) 

#do the same for the test cohort and validation cohorts separeatly

df.test = subset(df, Group=="Test")

ggplot(subset(df.test, Type != "HEALTHY"), aes(y=NE_meth_score, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_NEPC_test.pdf"), height=3, width=3)

ggplot(subset(df.test, Type != "HEALTHY"), aes(y=PRAD_meth_score, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_PRAD_test.pdf"), height=3, width=3)

ggplot(subset(df.test, Type != "HEALTHY"), aes(y=portion, x=Type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_ratio_test.pdf"), height=3, width=3)

# -----

message("test cohort")
p=wilcox.test(x=df.test$NE_meth_score[df.test$Type=="NEPC"], y=df.test$NE_meth_score[df.test$Type=="PRAD"])$p.value    
message("Mann Whitney test for NEPC vs PRAD NE_meth_scores: ", p) 

p=wilcox.test(x=df.test$PRAD_meth_score[df.test$Type=="NEPC"], y=df.test$PRAD_meth_score[df.test$Type=="PRAD"])$p.value                
message("Mann Whitney test for NEPC vs PRAD PRAD_meth_scores: ", p) 

p=wilcox.test(x=df.test$portion[df.test$Type=="NEPC"], y=df.test$portion[df.test$Type=="PRAD"])$p.value                
message("Mann Whitney test for NEPC vs PRAD portion scores: ", p) 

# do the same for the validation cohort

df.val = subset(df, Group=="Validation")

ggplot(subset(df.val, Type != "HEALTHY"), aes(y=NE_meth_score, x=Type)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_NEPC_val.pdf"), height=3, width=3)

ggplot(subset(df.val, Type != "HEALTHY"), aes(y=PRAD_meth_score, x=Type)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_PRAD_val.pdf"), height=3, width=3)

ggplot(subset(df.val, Type != "HEALTHY"), aes(y=portion, x=Type)) +
  geom_boxplot(outlier.shape = NA) +                                                                           
  geom_jitter(aes(color=Type), height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols)
ggsave(file.path(rms_dir, "rms_score_type_ratio_val.pdf"), height=3, width=3)

# -----

message("Validaiton cohort")                                                                                    
p=wilcox.test(x=df.val$NE_meth_score[df.val$Type=="NEPC"], y=df.val$NE_meth_score[df.val$Type=="PRAD"])$p.value
message("Mann Whitney test for NEPC vs PRAD NE_meth_scores: ", p)                                              

p=wilcox.test(x=df.val$PRAD_meth_score[df.val$Type=="NEPC"], y=df.val$PRAD_meth_score[df.val$Type=="PRAD"])$p.value          
 
message("Mann Whitney test for NEPC vs PRAD PRAD_meth_scores: ", p)                                            

p=wilcox.test(x=df.val$portion[df.val$Type=="NEPC"], y=df.val$portion[df.val$Type=="PRAD"])$p.value        
message("Mann Whitney test for NEPC vs PRAD portion scores: ", p)  

# compare ctDNA content in test and validation cohorts
#cols.group <- c("Test" = "purple", "Validation" = "darkorange")

ggplot(df.test, aes(y=ctDNA, x=Type, color=Type)) + 
  geom_boxplot(outlier.shape = NA, color="black") +
  geom_jitter(height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols) +
  ylim(0,1)
ggsave(file.path(rms_dir, "ctDNA_test.pdf"), height=3, width=3)

p=wilcox.test(x=df.test$ctDNA[df.test$Type=="NEPC"], y=df.test$ctDNA[df$Type=="PRAD"])$p.value
message("Mann Whitney test for ctDNA of test and validation cohorts: ", p)

ggplot(df.val, aes(y=ctDNA, x=Type, color=Type)) +                      
  geom_boxplot(outlier.shape = NA, color="black") +                                                        
  geom_jitter(height = 0) +                                                                                
  theme_classic() +                                                                                        
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  scale_color_manual(values = cols) +                                                                      
  ylim(0,1)
ggsave(file.path(rms_dir, "ctDNA_val.pdf"), height=3, width=3)

p=wilcox.test(x=df.val$ctDNA[df.val$Type=="NEPC"], y=df.val$ctDNA[df.val$Type=="PRAD"])$p.value
message("Mann Whitney test for ctDNA of test and validation cohorts: ", p)       


#scatterplot of PRAD vs NEPC scores:
ggplot(subset(df, Type != "HEALTHY"), aes(y=NE_meth_score, x=PRAD_meth_score, col=Type, size=ctDNA)) +
  geom_point() +
  scale_color_manual(values = cols) +
  theme_classic()
ggsave(file.path(rms_dir, "PRAD_vs_NEPC_score.pdf"), height=4, width=5)

# plot auc for NEPC score for calssifying presence of ANY NE features:
aucsamples = subset(df, Type=="PRAD" | Type=="NEPC")

curves <-  calculate_roc(D=as.numeric(aucsamples$Type=="NEPC"), M=aucsamples$NE_meth_score, ci = FALSE)
auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=aucsamples, aes(d = as.numeric(aucsamples$Type=="NEPC"), m = NE_meth_score), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_NEPC.pdf"),
  width = 3, height = 3)

message("classification AUC for NEPC score: ", calc_auc(auc_plot)$AUC)

curves <-  calculate_roc(D=as.numeric(aucsamples$Type=="NEPC"), M=aucsamples$portion, ci = FALSE)
auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=aucsamples, aes(d = as.numeric(aucsamples$Type=="NEPC"), m = portion), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_portion.pdf"),
  width = 3, height = 3)

message("classification AUC for NEPC/PRAD score ratio: ", calc_auc(auc_plot)$AUC)

#plot portion scores for test and validation cohorts separately
curves <-  calculate_roc(D=as.numeric(df.test$Type=="NEPC"), M=df.test$portion, ci = FALSE)
auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_roc(data=df.test, aes(d = as.numeric(df.test$Type=="NEPC"), m = portion), labels = FALSE, color="darkgray")

ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_portion_test.pdf"),
  width = 3, height = 3)

message("classification AUC for NEPC/PRAD score ratio for test cohort: ", calc_auc(auc_plot)$AUC)

curves <-  calculate_roc(D=as.numeric(df.val$Type=="NEPC"), M=df.val$portion, ci = FALSE)
auc_plot <- ggplot() +                                                                                     
  geom_line(data = curves, aes(x=FPF, y = TPF)) +                                                          
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +             
  geom_roc(data=df.val, aes(d = as.numeric(df.val$Type=="NEPC"), m = portion), labels = FALSE, color="darkgray")    
                                                                                                           
ggsave(auc_plot, file = file.path(rms_dir, "auc_curve_portion_val.pdf"),                                  
  width = 3, height = 3)

message("classification AUC for NEPC/PRAD score ratio for validation cohort: ", calc_auc(auc_plot)$AUC)


#-----


#print ctDNA stats
ct = df$ctDNA
message("cohort ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))
ct = df$ctDNA[df$Type=="PRAD"]
message("PRAD ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))
ct = df$ctDNA[df$Type=="NEPC"]
message("NEPC ctDNA stats: min=", min(ct), " max=", max(ct), " median=", median(ct))

#plot portion score vs ctDNA estimate for NEPC and PRAD:
ggplot(df, aes(y=portion, x=ctDNA, col=Type)) +  
  geom_point() + theme_classic() + scale_color_manual(values = cols) +  
  geom_smooth(method="lm", col="blue", linetype="longdash", data=subset(df, Type=="PRAD")) +      
  geom_smooth(method="lm", col="red", linetype="longdash", data=subset(df, Type=="NEPC")) +             
ggsave(file.path(rms_dir, "rms_score_portion_ctDNA_scatter.pdf"), height=3, width=5) 

cor=cor.test(x=df$portion[df$Type=="PRAD"], y=df$ctDNA[df$Type=="PRAD"])
message("correlation for PRAD sample portion score and ctDNA = ", cor$estimate, "(p = ", cor$p.value, ")")

cor=cor.test(x=df$portion[df$Type=="NEPC"], y=df$ctDNA[df$Type=="NEPC"])
message("correlation for NEPC sample portion score and ctDNA = ", cor$estimate, "(p = ", cor$p.value, ")")

# OS plots:
df$OS = as.numeric(df$OS)


#plot portion score vs ctDNA estimate for NEPC and PRAD:                                                         
ggplot(df, aes(y=portion, x=OS, col=Type)) +                                                       
  geom_point() + theme_classic() + scale_color_manual(values = cols) + 
  geom_smooth(method="lm", col="gray", linetype="longdash") 
ggsave(file.path(rms_dir, "rms_score_portion_OS.pdf"), height=3, width=5)                                         
                                     
cor=cor.test(x=df$portion[!is.na(df$OS)], y=df$OS[!is.na(df$OS)])      
message("correlation for portion score and OS = ", cor$estimate, "(p = ", cor$p.value, ")")                   

#same for ctDNA vs OS:
ggplot(df, aes(y=ctDNA, x=OS, col=Type)) +               
  geom_point() + theme_classic() + scale_color_manual(values = cols) + 
  geom_smooth(method="lm", col="gray", linetype="longdash")
ggsave(file.path(rms_dir, "ctDNA_OS.pdf"), height=3, width=5)                                         
                                                         
cor=cor.test(x=df$ctDNA[!is.na(df$OS)], y=df$OS[!is.na(df$OS)])                                             
message("correlation for ctDNA and OS = ", cor$estimate, "(p = ", cor$p.value, ")") 


#look at NEPC score as a function of ctDNA quantiles:
df$ctDNA_quantile = rep(NA, nrow(df))
cuts = quantile(df$ctDNA)

i=1
for(c in cuts) {
	df$ctDNA_quantile[df$ctDNA > c] = names(cuts[i])
	i=i+1
}
df$ctDNA_quantile = factor(df$ctDNA_quantile, levels = names(cuts))


ggplot(subset(df, Type != "HEALTHY" & !is.na(ctDNA_quantile)), aes(y=portion, col=ctDNA_quantile, x=Type)) +
  geom_boxplot(outlier.shape = NA) +                                                                           
  geom_jitter(color="gray", height = 0) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1)) +

ggsave(file.path(rms_dir, "rms_ctDNA.pdf"), height=3, width=6)


