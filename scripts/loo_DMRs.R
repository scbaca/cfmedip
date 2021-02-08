# Reads in MEDIPS objects, identifies DMRs, 
# and uses them to classify a held out sample for validation

# Based on script by Keegan Korthauer

# this is derived from DMRs.R, but uses a leave one out approach instead of train-test splits

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
source("scripts/limmamedips.R")

# filelist = string with comma-separated list of medip object files
# lab1 = label for cases
# lab2 = label for controls
# sig.level = qvalue threshold for defining DMRs (optional, not used when top is specified)
# mdip.opt = logical whether to adjust for gc content bias in counts
# chrs = vector of chromosome names to include
# regulatory = logical whether to restrict to regulatory regions
# blacklist = bed file specifying regions to remove from consideration 
# top = threshold for number of DMRs (take top `top` DMRs ranked by qval) - will be appended to results files (heatmaps, results tables)
# colnames = logical whether to include sample names in heatmaps
# out.dir = directory location to save results
# sample_name = name of held out sample to be used for validation
# iteration = a number that will be used to label subdirectory if running multiple iterations
# restrict_up = DMRs in the training data must overlap with this bed file and be UP in cases compared to controls
# restrict_down = DMRs in the training data must overlap with this bed file and be DOWN in cases compared to controls. Only used if restrict_up is set.
# depth_file = file with depth for each sample, will be included as a covariate in limma model if provided

# differential coverage
compute.diff <- function(filelist=files,
			metasheet=metasheet,
			lab1 = NULL, lab2 = NULL,
			sig.level = NULL,
			mdip.opt = FALSE,
			chrs =  paste0("chr", c(1:22, "X", "Y", "M")),
			regulatory = FALSE,
			blacklist = NULL,
			top = 300,
			colnames = TRUE,
			out.dir = "out", 
			sample_name = NULL,
			restrict_up = NULL,
			restrict_down = NULL,
			alpha = 0.1,
			depth_file = NULL){

	set.seed(1)

	if (!is.null(sample_name)) {
		out.dir = file.path(out.dir, "leave_one_out", sample_name)	
		dir.create(out.dir, recursive=T)
	}


        #read metasheet	
	met=read.table(metasheet, header=T, sep=",", as.is=T)
	if(!all(c("SampleName", "Class", "Type") %in% colnames(met))) stop("check metasheet file for required columns")

	#will batch info be included?
	if("Batch" %in% colnames(met)) show.batch=T

	#only use samples labeled case or control in metasheet
	met = subset(met, Class %in% c("case", "control"))

	#match files to samples in metasheet
        file=unlist(strsplit(filelist," "))
	SampleName = file %>% str_replace(".*/","") %>% str_replace(".medip.rds","") 
	if(!(all(met$SampleName %in% SampleName))) {
		message(met$SampleName[!met$SampleName %in% SampleName])
		stop("metasheet file does not match sample names")                     
	}
                                                                                                         
	message("excluding the following samples that are in the config file but not the metasheet file (if any):")       
	message(SampleName[!SampleName %in% met$SampleName])
	met = merge(met, data.frame(cbind(SampleName, file)))
	met$file = as.character(met$file)                                                                                                      
        # load medips objects:                                                                                           
        message("reading in medips files") 
        obj1 <- lapply(met$file[met$Class=="case"], readRDS) 
        obj2 <- lapply(met$file[met$Class=="control"], readRDS) 
	validate = FALSE

	# the coupling set contains counts of CG for each window, can be used for normalization
        if(file.exists("analysis/couplingset/couplingset.rds")) {
                CS = readRDS("analysis/couplingset/couplingset.rds")
        } else {
        	CS = MEDIPS.couplingVector(pattern = "CG", refObj = obj1[[1]]) 
        }

        # set graph labels and file names
	diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".diff.rds"))
	bed.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".diff.bed"))
	coef.bed.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".coef.bed"))
	model.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".model.rds"))
	plot.model.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".model.pdf"))
	plot.q.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".q.lfc.pdf"))
	dmrs.file <- file.path(out.dir, paste0(lab1, ".", lab2, ".dmrs.rds"))
	heatmap.file <- file.path(out.dir,paste0(lab1, ".", lab2, 
	  ".diff.heatmap.top", top,".pdf"))
	dmrs_val.file <- file.path(out.dir,
	  paste0(lab1, ".", lab2, ".dmrs_val.rds"))
	summary.file <- paste0(out.dir,"/performance.tsv")
	prob.table.file <- paste0(out.dir,"/sample_prob_table.tsv")
	diff.file <- file.path(out.dir, paste0(lab1, ".", lab2, 
	  ".diff.loo.rds"))
	heatmap.file <- file.path(out.dir, 
	  paste0(lab1, ".", lab2, ".diff.loo.heatmap.top", top,".pdf"))
	heatmap.file.test <- file.path(out.dir, 
	  paste0(lab1, ".", lab2, ".diff.loo.test.heatmap.top", top,".pdf"))

	#exclude the held out sample from training:

        # identify DMRs:
        message("===== Identifying DMRs in training set =====")

        message("Using the following cases:")
        l1 = unlist(lapply(obj1, function(x) {x@sample_name})) %>% str_replace(".dedup.bam", "") 
        message(paste(l1,delim=" "))                                                                                                  
        message("Using the following controls:")                                      
        l2 = unlist(lapply(obj2, function(x) {x@sample_name})) %>% str_replace(".dedup.bam", "") 
        message(paste(l2,delim=" ")) 


	ix1 <- !(l1 %in% sample_name) 
	ix2 <- !(l2 %in% sample_name) 

	n1=sum(ix1)
	n2=sum(ix2)

	if (!is.null(sig.level)){
		message("Using ", sig.level, " FDR cutoff instead of taking top ",
		  top, " regions.")
	}else{
		message("Using top ", top , " regions instead of FDR cutoff")
	}
	# uses a custom version of MEDIP.meth employing limma rather than ttest/edgeR

	if (!is.null(depth_file)){
		message("Including depth as covariates in the model")
		depth <- read.table(depth_file, stringsAsFactors=F, header=T)	
		#put depth dataframe in the same order as the medips objects
		depth <- depth[match(c(l1, l2), depth$sample),]
	} else {
		depth = NULL
	}
	diff = MEDIPS.meth(MSet1 = obj1[ix1], MSet2 = obj2[ix2],
	  CSet = CS, diff.method = "limma", chr = chrs,
	  p.adj = "BH", MeDIP = mdip.opt, minRowSum = 0.2*(n1+n2),
	  depth = depth)
	saveRDS(diff, file = diff.file)
	# diff has counts, rpkm, and adj p value for each window

	#NOTE: the way diff is calculated by MEDIPS.meth above, positive logFC means higher in CONTROLS (ie, MSet2 relative to MSet1)
	# reverse sign of logFC
	diff$logFC = -1 * diff$logFC

	# restrict to regulatory regions
	incl = !is.na(diff$P.Value)
	message("Trimming search regions. Starting with ", sum(incl), " windows")

	#restrict to windows with coupling factor > 0
	incl = incl & (CS@genome_CF > 0)
	message (sum(incl), " windows remain after restricting to windows with coupling factor > 0")

	if (regulatory){
      		# get cpg islands, shelves, and shores, plus enhancers
		annots <- c("hg19_cpgs", "hg19_enhancers_fantom")
		annotations = build_annotations(genome = 'hg19', annotations = annots)
		annotations <- annotations[annotations$type != "hg19_cpg_inter",]
		incl = incl & (makeGRangesFromDataFrame(diff) %over% annotations)
		message(sum(incl), " windows remain after restricting to CpGs/enhancers")
	}

	# restrict to specified sites (ie NEPC-up in LuCaPs):
	if (restrict_up != "none"){
		sites_up=import(restrict_up, format="BED") 
		incl_restrict = incl & (makeGRangesFromDataFrame(diff) %over% sites_up) & (diff$logFC > 0) 
		message(sum(incl_restrict), " windows remain after limiting to regions in restrict_up file")
		if (restrict_down != "none"){
			sites_down=import(restrict_down, format="BED")
                	incl_restrict = (incl_restrict | incl & (makeGRangesFromDataFrame(diff) %over% sites_down) & (diff$logFC < 0))
                	message(sum(incl_restrict), " windows remain after limiting to regions in restrict_down file")
		}
		incl = incl_restrict
	}
	

	if (!is.null(blacklist)){
		bl=import(blacklist, format="BED")
		incl = incl & !(makeGRangesFromDataFrame(diff) %over% bl)
		message(sum(incl), " windows remain after removing blacklisted regions")
	}

	diff = diff[incl,]

	# recalculate the adjusted p values for differential methylation (originally calculated for all windows, not restrictedd set)
	diff$limma.adj.p.value = p.adjust(diff$P.Value, "BH")

        # plot p val vs logFC	
	message("plotting pval vs logFC")
	plot = diff[!is.na(diff$limma.adj.p.value),]
	ggplot(data=plot[sample(1:nrow(plot),10000, replace=TRUE),])  + geom_point(aes(y=-1*log(limma.adj.p.value,10),x=logFC)) +
	  ylab("-log10 q-value") + xlab("log2 fold change (pos = higher in cases)") + theme_classic()
	ggsave(plot.q.file)

	# create a bed file with significant DMRs at 0.05 (these are not necessarily the final ones used for classification).
	message("saving bed file with DMRs")
	d=diff[diff$limma.adj.p.value < 0.05 & !is.na(diff$limma.adj.p.value),]
	bed=data.frame(chr=d$chr,start=d$start,
	  end=d$stop,case=d$MSets1.rpkm.mean,
	  control=d$MSets2.rpkm.mean, LFC=d$logFC,Padj=d$limma.adj.p.value)
	write.table(bed,file=bed.file,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

        #also create bed files with only up regions and only down regions 
	write.table(bed[bed$LFC>0,1:3], paste0(bed.file,".case.up"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE) 
	write.table(bed[bed$LFC<0,1:3], paste0(bed.file,".case.down"),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE) 
	colnames(diff)[3] <- "end" #change from default "stop"   

	# make heatmap of medip signal at dmrs
	if (!is.null(sig.level)){
		which.sig <- which(diff$limma.adj.p.value < sig.level & 
		               abs(diff$logFC) > 2 &
		               !is.na(diff$P.Value) )
	}else{# take top ntop

		which.sig <- which(rank(diff$limma.adj.p.value,
		  ties.method = "random") <= as.numeric(top))
	}
	which.up <- which(diff$logFC > 0)
	which.down <- which(diff$logFC < 0)
	
	coef.pos = which.sig %in% which.up
	coef.lower.lim = ifelse(coef.pos, 0, -Inf)
	coef.upper.lim = ifelse(coef.pos, Inf, 0)

	# get read counts for DMRs in each sample:
	dmrs <- diff[which.sig,
	  grepl("counts", colnames(diff))]
	dmrs <- dmrs[,!grepl("mean", colnames(dmrs))]

	rownames(dmrs) <- paste0(diff[which.sig,]$chr, ":", 
	  diff[which.sig,]$start, "-",
	  diff[which.sig,]$end)

	# adjust sample labels
	colnames(dmrs)[1:n1] <- paste0(lab1, "_", colnames(dmrs)[1:n1])
	colnames(dmrs)[(n1 + 1):(n1 + n2)] <- paste0(lab2, 
	  "_", colnames(dmrs)[(n1 + 1):(n1 + n2)])
	colnames(dmrs) <- gsub(".counts|.rpkm", "", colnames(dmrs))
	colnames(dmrs) <- gsub(".sorted.bam", "", colnames(dmrs))

	# normalize counts

	# select counts columns for training samples
	counts.train <- diff[, grepl("counts", colnames(diff))]
	counts.train <- counts.train[, !grepl("mean", colnames(counts.train))]

	# only use windows with >= 10 counts across samples for normalization
	which.norm <- which(rowSums(counts.train)>=10)

	d_train <- DGEList(counts=counts.train) # (makes an edgeR counts matrix (DGEList))

	# calculate normalization factor for each sample using TMM method, 
	# which removes the top and bottom x percentiles then normalizes to the mean 
	#(assumes that most regions are NOT differentially methylated)
	d_train <- calcNormFactors(d_train[which.norm,], refColumn = 1)

	# normalize (divide counts in each window by norm factor * library size / 1e6):
	dmrs <- sweep(dmrs, MARGIN=2, 
	  (d_train$samples$norm.factors*d_train$samples$lib.size)/1e6, `/`)

	saveRDS(dmrs, file = dmrs.file)

        #create heatmap of counts per million on log scale:
	theme_set(theme_bw())
	ecolors <- c("white",colorRampPalette( (brewer.pal(9, "Reds")) )(255)) # set colors
	type <- c("lightgrey", "black") #colors for heatmap labels
	names(type) <- c(lab1, lab2)
	ha_column = HeatmapAnnotation(df = data.frame(Type = ifelse(grepl(lab1, colnames(dmrs)), lab1, lab2)),
	   col = list(Type = type)) # heatmap column color annotations

	ht = Heatmap(log(dmrs+1), name = "log(CPM+1)", 
	  top_annotation = ha_column, col = ecolors,
	  show_row_names = FALSE, show_column_names = colnames,
	  column_names_gp = gpar(fontsize = 9),
	  column_title = paste0("Top ", top)) 

	pdf(heatmap.file, width=8)
	draw(ht)	
	dev.off()

	# validation
        message("===== Classifying held out sample =====")
	grp1 <- gsub(paste0(lab1, "_"), "", 
	  colnames(dmrs)[which(grepl(lab1, colnames(dmrs)))])
	grp2 <- gsub(paste0(lab2, "_"), "", 
	  colnames(dmrs)[which(grepl(lab2, colnames(dmrs)))])

	if (any(l1 %in% sample_name)) {
		m1 = MEDIPS.meth(MSet1 = obj1[which(l1 %in% sample_name)])[incl,]
	} else { 
		m1 = MEDIPS.meth(MSet1 = obj2[which(l2 %in% sample_name)])[incl,]
	}

	m1 <- m1[,!grepl("mean", colnames(m1)), drop = FALSE]
	counts.test <- m1[,grepl("counts", colnames(m1)), drop = FALSE]
	dmrs_new <- counts.test[which.sig,]  # counts for test set at DMRs defined by training set (not yet normalized)
	d_all <- DGEList(counts=cbind(counts.train, counts.test))
	d_all <- calcNormFactors(d_all[which.norm,], refColumn = 1)
	dmrs_new = dmrs_new / (d_all$samples$norm.factors[n1+n2+1] * d_all$samples$lib.size[n1+n2+1] / 1e6)

	#need to change dmrs_new back into a dataframe:
	dmrs_new = data.frame(dmrs_new)
	colnames(dmrs_new)=sample_name
	rownames(dmrs_new) <- paste0(diff[which.sig,]$chr, ":", 
	diff[which.sig,]$start, "-",
	diff[which.sig,]$end)

	saveRDS(dmrs_new, file = dmrs_val.file)

	rm(m1)
	
	# classify validation sample

	library(glmnet)
	#create a cross validation glmnet object 	
	### previously this was in a try-catch statement, put it back in one if needed
	cvob1 = cv.glmnet(x=t(log(dmrs+1)),y=c(rep(1, n1), rep(0, n2)),
	  family="binomial", alpha=alpha, type.measure="auc", 
	  nfolds = 3, lambda = 10^seq(2,-5,length.out = 100), 
	  upper.limits = coef.upper.lim, lower.limits = coef.lower.lim,
	  standardize=FALSE)

	# plot model parameters
	saveRDS(cvob1, file = model.file)
	pdf(plot.model.file)
	plot(cvob1$glmnet.fit,xvar="dev")
	plot(cvob1$glmnet.fit,xvar="lambda")
	plot(cvob1$glmnet.fit,xvar="norm")
	plot(cvob1)
	dev.off()	

	#save DMRs with their model coefficients to a bed file
	dmr_bed = data.frame(diff[which.sig,][c("chr", "start", "end", "MSets1.rpkm.mean", "MSets2.rpkm.mean", "logFC", "limma.adj.p.value")], 
	  coef = as.numeric(coef(cvob1)[2:nrow(coef(cvob1)),]))
	write.table(dmr_bed, file=coef.bed.file, quote = FALSE, sep = "\t", col.names=FALSE, row.names=FALSE)

	message("evaluating test sample")
	# assign test sample to groups using model generated on the training data above
	# "new" is the probability that is reported in the final output
	new <- predict(cvob1, newx=t(log(dmrs_new+1)), type = "response", 
	  s = "lambda.min") 

	#create df with true group labels, probabilities, and auc.
        ret_tab <- data.frame(sample_name = colnames(dmrs_new),                        
          true_label = ifelse(sample_name %in% l1, lab1, lab2),
          class_prob = as.vector(new),                                                 
          auc = NA) 

	mns <- data.frame(grp1 = rowMeans(log(dmrs[,(1:n1)]+1)),
	  grp2 = rowMeans(log(dmrs[,(n1+1):(n2+n1)]+1)))
	colnames(mns) <- c(lab1, lab2)
	newdist <- as.matrix(dist(t(cbind(log(dmrs_new+1), mns)))) # a distance matrix of each sample as well as means for each group (Euclidean by default)

	ret_tab$class_label <-  rep(NA, ncol(dmrs_new))

	# assign the sample to the closer group based on euclidean distance to the mean value for each group:
	ret_tab$class_label[1] <-                                            
	  ifelse(newdist[1,ncol(newdist)] < newdist[1,ncol(newdist)-1],    
	  lab2, lab1)  

	message("Held out sample ", colnames(dmrs_new)[1], 
	  " classified to group ", ret_tab$class_label[1], "(prob= ", 
	  ret_tab$class_prob[1], ")")

	message("Saving sample-level probability table")
	write.table(ret_tab, quote=FALSE, row.names=FALSE,
	  file= prob.table.file, 
	  sep = "\t")
}

args <- commandArgs( trailingOnly = TRUE )
filelist = args[1]
out.dir = args[2]
metasheet = args[3]
blacklist = args[4]
sample_name = args[5]
restrict_up = args[6]
restrict_down = args[7]
depth_file = args[8]

compute.diff(filelist=filelist, metasheet=metasheet,
  lab1 = "case", lab2 = "control", out.dir = out.dir,
  regulatory=FALSE, top=1000, blacklist = blacklist, sample_name = sample_name,
  alpha=0.1, restrict_up = restrict_up, restrict_down = restrict_down, depth_file=NULL)
