#create summary plots, including AUC plots and mean methylation "risk score" for each sample.

library(ggplot2)
library(dplyr)
library(readr)
library(plotROC)

args = commandArgs( trailingOnly = TRUE )

out.dir=args[1]
metasheet=args[2]
iter=args[3]
# --- tmp - will remove this once incorporated into snakemake workflow:
#out.dir="analysis/dmrs/cf_nepc_vs_prad"
#out.dir="analysis/dmrs/lucap_cf_nepc_vs_prad"
#out.dir="analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_hiconf"
#out.dir="analysis/dmrs/lucap_cf_nepc_vs_prad_train_test"
#out.dir="analysis/dmrs/lucap_cf_nepc_vs_prad_train_test_batch6"
#metasheet="run_files/metasheet.csv.lucap_cf_nepc_vs_prad_train_test_hiconf"
#metasheet="run_files/metasheet.csv.lucap_cf_nepc_vs_prad_train_test"
#metasheet="metasheet.csv"
# --- end tmp

# read sample annotations in metasheet file
meta=read.table(metasheet, header=T, sep=",")

files = c()
# for(i in sprintf("%03d", 1:100)){files <- c(files,paste0(out.dir,"/iter.",i,"/sample_prob_table.tsv"))}
# for(i in sprintf("%03d", 1:50)){files <- c(files,paste0(out.dir,"/iter.",i,"/sample_prob_table.tsv"))} #tmp
for(i in sprintf("%03d", 1:iter)){files <- c(files,paste0(out.dir,"/iter.",i,"/sample_prob_table.tsv"))} #tmp

theme_set(theme_bw())

tmp <- files %>%
  purrr::map(read_tsv)

tmp <- lapply(seq_along(files), 
function(x) {
    mutate(tmp[[x]], 
      filename=files[x],
      idx=paste0(true_label,1:nrow(tmp[[x]])))
  }) 

tmp <- tmp %>%
  do.call("rbind", .) 
idx <- stringr::str_locate(tmp$filename[1], "iter.")[2]
tmp$iteration = substr(tmp$filename, idx+1, idx+3)

tmp$SampleName = tmp$sample_name %>% 
  stringr::str_replace(".dedup.bam","") %>% 
#  stringr::str_replace(".dedup.bam.counts","") %>% 
  stringr::str_replace("case_","") %>% 
  stringr::str_replace("control_","")

tmp = merge(tmp, meta[,c("SampleName", "Type")])

#auc table
auc_summary <- tmp %>%
  group_by(iteration,Type) %>%
#  group_by(iteration,Stage) %>%
  summarize(auc=unique(auc)) %>%
  group_by(Type) %>%
#  group_by(Stage) %>%
  summarize(meanAUC = mean(auc, na.rm = TRUE),
            sdAUC = sd(auc, na.rm = TRUE),
            n = n(),
            lowerAUC = max(0, meanAUC - qnorm(0.975)*sdAUC/sqrt(n)),
            upperAUC = min(meanAUC + qnorm(0.975)*sdAUC/sqrt(n), 1))
auc_summary


#tmp
auc_summary <- tmp %>%                                        
  summarize(meanAUC = mean(auc, na.rm = TRUE),                
            sdAUC = sd(auc, na.rm = TRUE),                    
            n = n(),                                          
            lowerAUC = max(0, meanAUC - qnorm(0.975)*sdAUC/sqrt(n)),        
            upperAUC = min(meanAUC + qnorm(0.975)*sdAUC/sqrt(n), 1))        
auc_summary  
#end tmp

#(need to define out.dir)
write.table(format(data.frame(auc_summary), digits=3), quote=FALSE, row.names=FALSE, sep="\t",
  file=file.path(out.dir, "auc_summary_100iter.txt"))


#plot of auc
#cols <- c("met_UC" = "red", "control" = "blue", "localized_UC"="purple")
#cols <- c("metastatic" = "red", "MIBC" = "purple", "non-MIBC" = "blue", "control" = "green", "recurrent" = "orange")#, "non-recurrent" = "darkblue")
cols <- c("PRAD" = "blue", "NEPC" = "orange")
#cols <- c("Progression" = "red", "No Progression" = "blue") #tmp


tmp %>% 
#  group_by(iteration,Stage) %>%
  group_by(iteration,Type) %>%
  summarize(auc = mean(auc, na.rm = TRUE)) %>%
  ggplot(aes(x = Type, y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = cols) +
  geom_jitter(aes(color=Type), height = 0, alpha = 0.8)+
  ylim(0,1) +
  xlab("")+ylab("AUC") +
  theme(legend.position = "none")
ggsave(file=file.path(out.dir, "boxplot_auc_summary_100iter.pdf"),
  width=4, height=3)



# Reorder
#######################################
# Functions for sorting factor levels #
# (janhove.github.io)                 #
#######################################
# Sort factor levels by the factor level mean of another covariate
sortLvlsByVar.fnc <- function(oldFactor, sortingVariable, ascending = TRUE) {
  
  require("dplyr")
  require("magrittr")
  
  # Combine into data frame
  df <- data.frame(oldFactor, sortingVariable)
  
  ###
  ### If you want to sort the levels by, say, the median, sd etc. instead of the mean,
  ### just change 'mean(sortingVariable)' below to, say, 'median(sortingVariable)'.
  ###
  
  # Compute average of sortingVariable and arrange (ascending)
  if (ascending == TRUE) {
    df_av <- df %>% group_by(oldFactor) %>% 
      summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(meanSortingVariable)
  }
  
  # Compute average of sortingVariable and arrange (descending)
  if (ascending == FALSE) {
    df_av <- df %>% group_by(oldFactor) %>%
      summarise(meanSortingVariable = median(sortingVariable, na.rm=TRUE)) %>% 
      arrange(desc(meanSortingVariable))
  }
  
  # Return factor with new level order
  newFactor <- factor(oldFactor, levels = df_av$oldFactor)
  return(newFactor)
}
####  end function to reorder


risk_summary <- tmp %>%
  group_by(sample_name, Type) %>%
#  group_by(sample_name, Stage) %>%
  summarize(mean_risk_score = mean(class_prob),
            sd_risk_score = sd(class_prob, na.rm = TRUE),
            n = n(),
            lower_risk_score = mean_risk_score - qnorm(0.975)*sd_risk_score/sqrt(n),
            upper_risk_score = mean_risk_score + qnorm(0.975)*sd_risk_score/sqrt(n),
            true_label = unique(true_label))

write.table(format(data.frame(risk_summary), digits=2), quote=FALSE, row.names=FALSE, sep="\t",
  file=file.path(out.dir, "risk_score_summary.txt"))


sample_probs_all <- tmp 

sub <- sample_probs_all %>% mutate(truth = as.numeric(true_label=="case"))
curves <- sub %>%
  group_by(iteration) %>%
  group_map(~ calculate_roc(D=.x$truth, M=.x$class_prob, ci = FALSE, alpha = 0.05))
curves <- lapply(1:length(files), function(x) {
  cbind(curves[[x]], iteration=x)
})
curves <- do.call(rbind, curves)



ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF, group=iteration ), alpha=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  geom_roc(data=sub, aes(d = truth, m = class_prob), labels = FALSE, color="blue") 
ggsave(file =file.path(out.dir, "auc_curves.pdf"),
  width = 3, height = 3)



sample_probs <- sample_probs_all
sample_probs$sample_name <- as.factor(sample_probs$sample_name)
sample_probs$sample_name <- sortLvlsByVar.fnc(sample_probs$sample_name,
  sample_probs$class_prob)


p <- sample_probs %>%
  ggplot(aes(x = sample_name, y = class_prob, color = Type)) +
#  ggplot(aes(x = sample_name, y = class_prob, color = Stage)) +
  geom_boxplot() +
#  theme(axis.text.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #tmp
  labs(color = "Sample Type") + 
  scale_color_manual(values = cols) + 
  ylab("Methylation score") + xlab("Samples") +
  ylim(0,1)
ggsave(p,
  file =file.path(out.dir, "boxplot_risk_score_summary.pdf"),
  width = 12, height = 4)
#  width = 6, height = 3)


