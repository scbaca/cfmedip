# This script creates MEDIPS objects and carries out various QC and 
# exploratory plotting procedures

# Based on script by Keegan Korthauer

# call as: "Rscript scripts/MeDIPS.R {sample} {bam} {output directory}

suppressMessages(library(MEDIPS))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
source("scripts/MEDIPS.CpGenrich.new.R") #an updated MEDIPS.CpGenrich function that doesn't require RangedData, which is no longer supported
medips_function <- function(sample, bam, outdir, genome.relH.file, genome.GoGe.file) {

        set.seed(1)

	# medip object parameters	
	BSgenome="BSgenome.Hsapiens.UCSC.hg19"
	uniq=1e-3 # p-value for determining max number of reads at some position to be allowed
	extend=300 # extend reads by this length (ie, to reach the estimated insert size)
	shift=0 # shift reads by this length
	ws=300 #window size in base-pairs for binning MeDIP-seq reads. May make this a parameter in the config file eventually

	# canonical chrs
	chr.select <- paste0("chr", c(1:22, "X", "Y", "M"))

        # CpG enrichment - this runs slowly. In future may want to get this data manually from fastqs rather than using MEDIPS to get it from bams.
       er <- MEDIPS.CpGenrich.new(file = bam, BSgenome = BSgenome,
	  chr.select = chr.select, extend = extend, shift = shift, 
	  uniq = uniq, paired = F, oneoutof=5000, 
	  genome.relH.file = genome.relH.file, genome.GoGe.file = genome.GoGe.file)
	sink(file.path(outdir, paste0(sample, ".CpG.enrichment.txt")))
	cat(sample, "enrichment.score.relH:", er$enrichment.score.relH,"\n") 
        cat(sample, "enrichment.score.GoGe:", er$enrichment.score.GoGe,"\n") #will add a script to aggregate these data across samples
	sink()

	# create medip objects
	medip = MEDIPS.createSet(file = bam, BSgenome = BSgenome,
	 extend = extend, shift = shift, uniq = uniq,
	 window_size = ws, chr.select = chr.select, paired = F)	  
	saveRDS(medip, file = file.path(outdir, paste0(sample, ".medip.rds")))

	# export wig file for medip object
	MEDIPS.exportWIG(Set=medip, file=file.path(outdir, paste0(sample,".wig")), format="rpkm", descr=sample)

	# "coupling set" with CG content of each window
	if(file.exists("analysis/couplingset/couplingset.rds")) {	
		CS = readRDS("analysis/couplingset/couplingset.rds")
	} else {
		CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip)
	}

	# saturation plots
	pdf(file.path(outdir, paste0(sample,".saturation.pdf"))) 

	sr = MEDIPS.saturation(file = bam, BSgenome = BSgenome,
	 uniq = uniq, extend = extend, shift = shift, window_size = ws,
	 chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
	 rank = FALSE)
	MEDIPS.plotSaturation(sr, main = paste0(sample, " Saturation analysis"))
	dev.off()
	
	# coverage plots
	pdf(file.path(outdir, paste0(sample,".CpGcoverage.pdf"))) 

	cr = MEDIPS.seqCoverage(file = bam, pattern = "CG",     
	 BSgenome = BSgenome, chr.select = chr.select, extend = extend,     
	 shift = shift, uniq = uniq)
	MEDIPS.plotSeqCoverage(seqCoverageObj=cr, type="pie", 
	 cov.level = c(0,1, 2, 3, 4, 5), main = paste0(sample, " CpG Coverage"))
	dev.off()

}

args <- commandArgs( trailingOnly = TRUE)
arg_sample = args[1]
arg_bam = args[2]
arg_outdir = args[3]
genome.relH.file=args[4]
genome.GoGe.file=args[5]

medips_function(arg_sample, arg_bam, arg_outdir, genome.relH.file, genome.GoGe.file)
