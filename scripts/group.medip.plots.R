# create various qc and analysis plots from medips objects
# Based on script by Keegan Korthauer

suppressMessages(library(MEDIPS))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(annotatr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggdendro))
suppressMessages(library(matrixStats))
suppressMessages(library(stringr))
suppressMessages(library(circlize))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(stringr))
# TODO: review whether all these libraries are still needed


group.plots <- function(filelist = NULL, metasheet = NULL,
			outdir = "out", regulatory = TRUE, 
			blacklist = blacklist, 
			restrict_up = NULL, 
			restrict_down = NULL,
			filter.ctDNA = TRUE) {

	set.seed(1)

	#read metasheet
	met=read.table(metasheet, header=T, sep=",", as.is=T)
	if(!all(c("SampleName", "Type") %in% colnames(met))) stop("check metasheet file for required columns")

	#will batch info be included?
#	if("Batch" %in% colnames(met)) show.batch=T
	show.batch=F

	#only use samples labeled case or control in metasheet. remove healthy samples
	met = subset(met, Class %in% c("case", "control") & Type != "HEALTHY")


	#filter by ctDNA content?
	if(filter.ctDNA) {
		met = subset(met, ctDNA >= 0.03) 
		message("filtering samples wtih ctDNA < 0.03")
	}

	#match files to samples in metasheet
	file=unlist(strsplit(filelist," ")) 

	SampleName = file %>% str_replace(".*/","") %>% str_replace(".medip.rds","")
	if(!(all(met$SampleName %in% SampleName))) stop("metasheet file does not match sample names")
	
	message("excluding the following samples that are in the config file but not the metasheet file (if any):")
	message(SampleName[!SampleName %in% met$SampleName])

	met = merge(met, data.frame(cbind(SampleName, file)))
	met$file = as.character(met$file)

	message("reading in medips files")
	medip = lapply(met$file, readRDS)

        if(file.exists("analysis/couplingset/couplingset.rds")) {
                CS = readRDS("analysis/couplingset/couplingset.rds")
        } else {
		message("generating coupling set")
        	CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip[[1]]) 
        }

	get_depths <- function(mdobjlist, samples){
	  depth <- data.frame(sapply(mdobjlist, function(x){
	  x@genome_count
	  }))
	  colnames(depth) <- samples
	  return(depth)
	}

	message ("getting read counts per window")
	df <- get_depths(medip, met$SampleName)

        grp = met$Type
	ids = met$SampleName

        # plot sparsity 
        message("creating sparsity plot")
	sparsity = colSums2(df == 0) / nrow(df)
	depth = colSums2(as.matrix(df))
	mtab <- data.frame(sample = met$SampleName, sparsity = sparsity, 
	  depth = depth,group = met$Type)                                         
	if(show.batch) mtab$batch = as.factor(met$Batch)                                              

	theme_set(theme_classic())
	if(show.batch) {
        mtab %>% ggplot(aes(x = sparsity, y = depth, color = group, shape = batch)) +                                                                
          geom_point() 
	} else {                                                                                                               
        mtab %>% ggplot(aes(x = sparsity, y = depth, color = group)) +                                                                
          geom_point() 
	}
        ggsave(file.path(outdir, "sparsity_depth.pdf"), height=4, width = 5)       

	sparsity_depth_file=file.path(outdir, "sparsity_depth.txt")
	message("writing sparsity and depth to ", sparsity_depth_file)
	write.table(mtab, file=sparsity_depth_file, row.names=F, col.names=T, quote=F, sep="\t")

	message("normalizing read counts")
	sf <- DESeq2::estimateSizeFactorsForMatrix(df)
	df <- sweep(df, MARGIN=2, sf, `/`)

	# restrict sites for pca/clustering as specified
	message("restricting sites as specified for pca and clustering")

	#get window coordinates to allow filtering based on bed files
        win=MEDIPS.meth(MSet1=medip[[1]], MeDIP=F) %>%
          select(c("chr", "start", "stop")) %>%
          makeGRangesFromDataFrame()
        incl = rep(TRUE, length(win))

	# restrict to regulatory regions
	message("Trimming search regions. Starting with ", sum(incl), " windows")

	#restrict to windows with coupling factor > 0
	incl = incl & (CS@genome_CF > 0)
	message (sum(incl), " windows remain after restricting to windows with coupling factor > 0")

	if (regulatory){
      		# get cpg islands, shelves, and shores, plus enhancers
		annots <- c("hg19_cpgs", "hg19_enhancers_fantom")
		annotations = build_annotations(genome = 'hg19', annotations = annots)
		annotations <- annotations[annotations$type != "hg19_cpg_inter",]
		incl = incl & (win %over% annotations)
		message(sum(incl), " windows remain after restricting to CpGs/enhancers")
	}

	# restrict to specified sites (ie NEPC-up in LuCaPs):
	if (restrict_up != "none" && restrict_up!= ""){
		sites_up=import(restrict_up, format="BED") 
		incl_restrict = incl & (win %over% sites_up)
		message(sum(incl_restrict), " windows remain after limiting to regions in restrict_up file")
		if (restrict_down!="none" && restrict_down!=""){
			sites_down=import(restrict_down, format="BED")
                	incl_restrict = (incl_restrict | incl & (win %over% sites_down))
                	message(sum(incl_restrict), " windows remain after limiting to regions in restrict_down file")
		}
		incl = incl_restrict
	}
	
	if (blacklist != "none"){
		bl=import(blacklist, format="BED")
		incl = incl & !(win %over% bl)
		message(sum(incl), " windows remain after removing blacklisted regions")
	}

	df = df[incl,]
	# require at least 20% of samples to have non-zero read counts for window to be included
	df <- df %>%
	filter(rowMeans(df > 0) > 0.20) 
	df <- df %>%
	filter(rowMeans(df) > 1)

	# pca plots
	message("generating pca plots")
	pcs <- prcomp(t(log(df+1)), center=TRUE, scale.=TRUE)

	#variance explained:
	eig = pcs$sdev^2
	varexp1 = 100*signif(eig[1]/sum(eig), digits=2)
	varexp2 = 100*signif(eig[2]/sum(eig), digits=2)
	varexp3 = 100*signif(eig[3]/sum(eig), digits=2)
	varexp4 = 100*signif(eig[4]/sum(eig), digits=2)

	tidydf <- select(data.frame(pcs$x), "PC1", "PC2", "PC3", "PC4") %>%
	  mutate(type = grp, id = ids)
	if(show.batch) tidydf$batch = as.factor(met$Batch)
	saveRDS(pcs,file.path(outdir,"pca.data.rds"))
	write.table(pcs$x, file=file.path(outdir,"pca.data.tsv"), row.names=T, col.names=T, quote=F, sep="\t")

	if(show.batch) {
	ggplot(tidydf, aes(x=PC1, y=PC2, colour = type, shape = batch)) + 
	  geom_point(size = 1) + 
	  xlab(paste0("PC1 (", varexp1, "% of variance explained)")) +
	  ylab(paste0("PC2 (", varexp2, "% of variance explained)"))
	ggsave(file.path(outdir, "PCA_1_2_cf_plasma.pdf"), width = 4, height = 3)

	ggplot(tidydf, aes(x=PC2, y=PC3, colour = type, shape = batch)) + 
	  geom_point(size = 1) + 
	  xlab(paste0("PC2 (", varexp2, "% of variance explained)")) +
	  ylab(paste0("PC3 (", varexp3, "% of variance explained)"))
	ggsave(file.path(outdir, "PCA_2_3_cf_plasma.pdf"), width = 4, height = 3)

	ggplot(tidydf, aes(x=PC3, y=PC4, colour = type, shape = batch)) +                              
	  geom_point(size = 1) +                                                        
	  xlab(paste0("PC3 (", varexp3, "% of variance explained)")) +
	  ylab(paste0("PC4 (", varexp4, "% of variance explained)"))
	ggsave(file.path(outdir, "PCA_3_4_cf_plasma.pdf"), width = 4, height = 3) 
	} else {
        ggplot(tidydf, aes(x=PC1, y=PC2, colour = type)) +
          geom_point(size = 1) +
	  xlab(paste0("PC1 (", varexp1, "% of variance explained)")) +
	  ylab(paste0("PC2 (", varexp2, "% of variance explained)"))
        ggsave(file.path(outdir, "PCA_1_2_cf_plasma.pdf"), width = 4, height = 3)

        ggplot(tidydf, aes(x=PC2, y=PC3, colour = type)) +
          geom_point(size = 1) +
	  xlab(paste0("PC2 (", varexp2, "% of variance explained)")) +
	  ylab(paste0("PC3 (", varexp3, "% of variance explained)"))
        ggsave(file.path(outdir, "PCA_2_3_cf_plasma.pdf"), width = 4, height = 3)                                                                  
        ggplot(tidydf, aes(x=PC3, y=PC4, colour = type)) + 
	   geom_point(size = 1) +
	   xlab(paste0("PC3 (", varexp3, "% of variance explained)")) +
	   ylab(paste0("PC4 (", varexp4, "% of variance explained)"))
        ggsave(file.path(outdir, "PCA_3_4_cf_plasma.pdf"), width = 4, height = 3) 
	}

	# TODO: add sample-sample correlation heatmap

	# dendrogram
	message("creating clustering plot")
	dm <- dist(t(log(as.matrix(df) + 1)))
	ggdendrogram(hclust(dm))
	hc <- hclust(dm)
	dendr <- dendro_data(hc, type="rectangle")
	clust.gr<-data.frame(label = colnames(df), type=grp)
	text.df<-merge(label(dendr),clust.gr,by="label")

	ggplot() +
	  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
	  geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=type),size=2) +
	  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
	  theme(axis.line.y=element_blank(),
	    axis.ticks.y=element_blank(),
	    axis.text.y=element_blank(),
	    axis.title.y=element_blank(),
	    panel.background=element_rect(fill="white"),
	    panel.grid=element_blank())
	ggsave(file.path(outdir, "hclust.pdf"), width = 7.5, height = 6.5)

	if(show.batch){
	text.df$batch = factor(met$Batch)
        ggplot() + 
          geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
          geom_text(data=text.df, aes(x=x, y=y, label=label, hjust=0,color=batch),size=2) + 
          coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
          theme(axis.line.y=element_blank(), 
            axis.ticks.y=element_blank(), 
            axis.text.y=element_blank(), 
            axis.title.y=element_blank(), 
            panel.background=element_rect(fill="white"), 
            panel.grid=element_blank()) 
        ggsave(file.path(outdir, "hclust.batch.pdf"), width = 7.5, height = 6.5)  
	}

	# how much coverage in windows with no CpGs?
	message("calculating coverage in windows with no CpGs")
	noCpGcov <- function(mdobjlist, CS){
		pctNoCpG <- sapply(mdobjlist, function(x){
			sum(x@genome_count[CS@genome_CF == 0]) / 
			  sum(x@genome_count)
  		})
		pctNonzeroCov <- sapply(mdobjlist, function(x){
  	  		sum(x@genome_count > 0) / length(x@genome_count)
  		})
		data.frame(pctNoCpG = pctNoCpG,
		 pctNonzeroCov = pctNonzeroCov)
	}

	#TODO: include this as an else in the if statement below
	df <- noCpGcov(medip, CS=CS)
	df$type = met$Type
	ggplot(df, aes(x = pctNoCpG*100)) +
          geom_histogram(breaks = seq(0,10,1)) +
	  xlim(0,10) +
	  facet_wrap(~type) +
	  xlab("Percent of total reads mapping to bins with 0 CpGs")
	ggsave(file.path(outdir, "PctReadsinBinswithZeroCpGs.pdf"))

	ggplot(df, aes(x = pctNonzeroCov*100)) +
	  geom_histogram(breaks = seq(0,100,10)) + 
	  xlim(0,100) +
	  xlab("Percent of bins with nonzero coverage")+
	  facet_wrap(~type) 
	ggsave(file.path(outdir, "PctBinswithNonzeroCov.pdf"))

	if(show.batch) {
 	df$batch = as.factor(met$Batch)
        ggplot(df, aes(x = pctNoCpG*100)) +
          geom_histogram(breaks = seq(0,10,1)) +
          xlim(0,10) +
          facet_wrap(~type + batch) +
          xlab("Percent of total reads mapping to bins with 0 CpGs")
        ggsave(file.path(outdir, "PctReadsinBinswithZeroCpGs.pdf"))

        ggplot(df, aes(x = pctNonzeroCov*100)) +
          geom_histogram(breaks = seq(0,100,10)) +
          xlim(0,100) +
          xlab("Percent of bins with nonzero coverage")+
          facet_wrap(~type + batch)
        ggsave(file.path(outdir, "PctBinswithNonzeroCov.pdf"))
	}

	# explore sample depths
	message("creating coverage plots")

	depths <- function(mdobjlist, CS, samples) {
		depth <- data.frame(sapply(mdobjlist, function(x){
			  x@genome_count[CS@genome_CF > 0]
       		}))
	  	colnames(depth) = samples
		depth <- depth %>% 
		filter(rowSums(depth == 0) < ncol(depth)) %>%
		tidyr::gather("sample", "count") 
		depth %>% filter(count > 0) %>%
	
		ggplot(aes(x = count)) +
		  geom_histogram(bins=10) +
		  facet_wrap(~sample) +
		  scale_x_log10() +
		  scale_y_log10() +
		  xlab("Coverage") +
		  theme(axis.text.x = element_text(angle = 90, hjust = 1))
		ggsave(file.path(outdir, paste0("NonzeroCovDistBySample.pdf")), 
		height = 9, width = 9)
        }
	depths(medip, CS, met$SampleName)
}

args <- commandArgs( trailingOnly = TRUE )
medip.files = args[1]
out.dir = args[2]
metasheet = args[3]
blacklist = args[4]
restrict_up = args[5]
restrict_down = args[6]


group.plots(filelist=medip.files, metasheet=metasheet,
  outdir = out.dir, blacklist = blacklist, 
  restrict_up = restrict_up, restrict_down = restrict_down, regulatory = FALSE)
