#plot methylation "risk score" for each sample. Based on plot.performance.R

library(ggplot2)
library(dplyr)
library(readr)
library(plotROC)

args = commandArgs( trailingOnly = TRUE )
out.dir=args[1]
metasheet=args[2]

# read sample annotations in metasheet file
meta=read.table(metasheet, header=T, sep=",")

# keep only samples labeled case or control
meta = subset(meta, Class %in% c("case", "control"))

files = c()
for (s in meta$SampleName){files <- c(files,file.path(out.dir, s, "sample_prob_table.tsv"))}

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

tmp$SampleName = tmp$sample_name %>% 
  stringr::str_replace(".dedup.bam","") %>% 
  stringr::str_replace("case_","") %>% 
  stringr::str_replace("control_","") %>%
  stringr::str_replace("\\.","-")
#the last line above is necessary because the '-' in the labels are replaced with "."

tmp = merge(tmp, meta[,c("SampleName", "Type")])

#plot of auc
cols <- c("PRAD" = "blue", "NEPC" = "orange")

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

sample_probs_all <- tmp 

sub <- sample_probs_all %>% mutate(truth = as.numeric(true_label=="case"))
curves <-  calculate_roc(D=sub$truth, M=sub$class_prob, ci = FALSE, alpha = 0.05)

auc_plot <- ggplot() +
  geom_line(data = curves, aes(x=FPF, y = TPF), alpha=0.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  geom_roc(data=sub, aes(d = truth, m = class_prob), labels = FALSE, color="blue")

auc_plot <- auc_plot + geom_text(label=as.character(calc_auc(auc_plot)$AUC), x=0.75, y=0.25) ### TODO: fix this, label doesn't display currently
ggsave(auc_plot, file =file.path(out.dir, "auc_curves.pdf"),
  width = 3, height = 3)

message("classification AUC: ", calc_auc(auc_plot)$AUC)

sample_probs <- sample_probs_all
sample_probs$sample_name <- as.factor(sample_probs$sample_name)
sample_probs$sample_name <- sortLvlsByVar.fnc(sample_probs$sample_name,
  sample_probs$class_prob)


p <- sample_probs %>%
  ggplot(aes(x = sample_name, y = class_prob, color = Type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(color = "Sample Type") + 
  scale_color_manual(values = cols) + 
  ylab("Methylation score") + xlab("Samples") +
  ylim(0,1)
ggsave(p,
  file =file.path(out.dir, "boxplot_risk_score_summary.pdf"),
  width = 12, height = 4)

#boxplots of PRAD and NEPC risk scores
p <- sample_probs %>%
  ggplot(aes(x = Type, y = class_prob)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Type), height = 0) +
  scale_color_manual(values = cols) +
  ylab("Methylation score") + xlab("Type") +
  ylim(0,1) +
  theme_classic()
ggsave(p,
  file =file.path(out.dir, "boxplot_risk_score_by_group.pdf"),
  width = 3, height = 3)

message("Mann Whitney test for NEPC risk scores by group:")
wilcox.test(x=sample_probs$class_prob[sample_probs$Type=="NEPC"], y=sample_probs$class_prob[sample_probs$Type=="PRAD"])
