##########################################################################
##Function calculates density of C's, G's, and CpGs in the given sequences as well as in the reference genome
##Returns CpG enrichment values for the given sequences.
##########################################################################
##Input:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##Param:	BSgenome, pattern (character), chromosomes, shift, extend
##Output:	CpGenrich results object
##Requires:	BSgenome, Biostrings, GenomicRanges, MEDIPS.Bam2GRanges, MEDIPS.txt2Granges, MEDIPS.GenomicCoordinates
##Modified:	11/10/2012
##Author:	Joern Dietrich, Lukas Chavez

#Updated 9/20/19 by SB to no longer require RangedData / IRanges::lapply, which are no longer supported
#Also now just analyzes a sample of reads to increase speed (set by "oneoutof")
#Can now pass rds files with the values of genome.relH and genome.GoGe so they don't need to be recalculated each time.

library(gtools)

MEDIPS.CpGenrich.new <-function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F, oneoutof=100, genome.relH.file="", genome.GoGe.file=""){

	## Proof correctness....
	if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
	
	## Read region file		
	fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
	path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/") 
	if(path==""){path=getwd()}		
	if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}	

	dataset = get(ls(paste("package:", BSgenome, sep = "")))	

	if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
	else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
	
	## Sort chromosomes
	if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
	if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}
	
	## Get chromosome lengths for all chromosomes within data set.
	cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))		

	chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])
	ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)
	
	##Calculate CpG density for regions
	total=length(chromosomes)
	cat("Calculating CpG density for given regions...\n")  
#	seq=matrix(unlist(IRanges::lapply(RangedData(GRange.Reads),function(x){		
#		i=which(mixedsort(chromosomes)%in%names(x) )
#		ranges(x)<-restrict(ranges(x),end=chr_lengths[which(chromosomes %in% names(x))])
#		y=DNAStringSet(getSeq(dataset, names=space(x), start=start(x), end=end(x), as.character=TRUE))
#		c(sum(as.numeric(vcountPattern("CG",y))),sum(as.numeric(vcountPattern("C",y))),sum(as.numeric(vcountPattern("G",y))),sum(as.numeric(width(y))),length(y))
#		}
#	),use.names=F),ncol=5,nrow=total,byrow=T)

#NEW 
#seq=matrix(unlist(lapply(seq_along(GRange.Reads),function(r){
nreads=ceiling(length(GRange.Reads)/oneoutof)
cat("sampling", nreads, "out of", length(GRange.Reads), "reads\n")
test=GRange.Reads #i=which(mixedsort(chromosomes)%in%seqnames(test[r]));
seq=matrix(unlist(lapply(seq(from=1, to=length(test), by=oneoutof),function(r){
#	test=GRange.Reads; #i=which(mixedsort(chromosomes)%in%seqnames(test[r]));
	ranges(test[r])<-restrict(ranges(test[r]),
	  end=chr_lengths[which(chromosomes %in% seqnames(test[r]))]); 
	y=DNAStringSet(getSeq(dataset,names=seqnames(test[r]), 
	  start=start(test[r]), end=end(test[r])));
	 # start=start(test[r]), end=end(test[r]), as.character=TRUE));
	c(as.numeric(vcountPattern("CG",y)),
	  as.numeric(vcountPattern("C",y)),
	  as.numeric(vcountPattern("G",y)),
	  width(y),length(y))
	}),use.names=F),ncol=5,nrow=nreads,byrow=T)
#        c(sum(as.numeric(vcountPattern("CG",y))),
#          sum(as.numeric(vcountPattern("C",y))),
#          sum(as.numeric(vcountPattern("G",y))),
#          sum(as.numeric(width(y))),length(y))
#        }),use.names=F),ncol=5,nrow=total,byrow=T)
	cat("done processing reads for CpG enrichment calculation\n")
	Value=colSums(seq)
#	unused=length(GRange.Reads)-Value[5]
#	if ( unused!=0 )cat(unused,"unused sequences, limits out of range\n")
	
	regions.CG=Value[1]
	regions.C=Value[2]
	regions.G=Value[3]
	all.genomic=Value[4]
	
	regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100	#percent of sequenced nucleotides that are CG
	regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))  
	

	##Calculate CpG density for reference genome
	if (genome.relH.file=="" || genome.GoGe.file=="") {

		cat(paste("Calculating CpG density for the reference genome", BSgenome, "...\n", sep = " "))	
		CG <- DNAStringSet("CG")
		pdict0 <- PDict(CG)
		params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
		genome.CG=sum(bsapply(params, pdict = pdict0))			
		params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
		alphabet=bsapply(params)
		genome.l=sum(as.numeric(alphabet))
 	
		genome.C=as.numeric(sum(alphabet[2,]))
		genome.G=as.numeric(sum(alphabet[3,]))
		genome.relH=genome.CG/genome.l*100
		genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);

	} else {
		cat(paste("Using the following files for CpG density calculations for the reference genome", genome.relH.file, genome.GoGe.file, "...\n", sep = " "))
		genome.relH=readRDS(genome.relH.file)
		genome.GoGe=readRDS(genome.GoGe.file)
	}

# modificaiton - save these values so they don't need to be calculated each time:
#	saveRDS(genome.relH,"genome.relH.rds")
#	saveRDS(genome.GoGe,"genome.GoGe.rds")
	
	enrichment.score.relH=regions.relH/genome.relH	
	enrichment.score.GoGe=regions.GoGe/genome.GoGe	
 	
	gc()
#	return(list(genome=BSgenome, regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G, regions.relH=regions.relH, regions.GoGe=regions.GoGe, genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG, genome.relH=genome.relH, genome.GoGe=genome.GoGe, enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))
	return(list(enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))
}
