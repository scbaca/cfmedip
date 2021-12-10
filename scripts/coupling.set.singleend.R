# Create a coupling set from the first sample for use by other rules
#call as: "Rscript scripts/coupling.set.R {input} {output}"

suppressMessages(library(MEDIPS))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

# may make these parameters in the config file eventually:
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1e-3
extend=300
shift=0
ws=300 
chr.select <- paste0("chr", c(1:22, "X", "Y", "M"))

coupling.set <- function(sample, output) {

	if(endsWith(sample, ".bam")) {       
		medip = MEDIPS.createSet(file = sample, BSgenome = BSgenome,
		  extend = extend, shift = shift, uniq = uniq,
		  window_size = ws, chr.select = chr.select, paired = F)
	} else { #assuming a medips input:
		medip = readRDS(sample)
	}
	CS = MEDIPS.couplingVector(pattern = "CG", refObj = medip)
	saveRDS(CS, file = output)
}

args <- commandArgs( trailingOnly = TRUE)
sample = args[1]
output = args[2]
coupling.set(sample=sample, output=output)
