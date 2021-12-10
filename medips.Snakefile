configfile: "config.yaml"

import pandas as pd
import yaml
import os

#load metasheet file containing sample metadata
metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')

# function to determine if input is a fastq file (as opposed to .medips file)
def input_is_fq(sample):
	first_file = config["samples"][sample][0]
	if first_file.endswith('.fq.gz'):
		return True
	else:
		return False

# get fq files with read 1 for paired-end reads
def get_fq1(wildcards):
	return config["samples"][wildcards.sample][0]

# get fq files with read 2 for paired-end reads
def get_fq2(wildcards):
	if(len(config["samples"][wildcards.sample])) == 1:
		return config["samples"][wildcards.sample][0]
		print("only one file provided for this sample, assuming it is a medips object")
	else:
		return config["samples"][wildcards.sample][1]

# --- functions to create lists of snakemake target files ---

# create targets for "all" rule
def get_targets(wildcards):
	ls = []
	
	#multiqc target:
	ls.append("analysis/multiqc/multiqc.html")

	#plot_qc target:
	ls.append("analysis/qc_plots/" + config['run'] + "/main_chr_unique_reads.pdf")
	ls.append("analysis/qc_plots/" + config['run'] + "/CpG_enrichment.pdf")

	#group_plots targets:
	ls.append("analysis/dmrs/" + config['run'] + "/hclust.pdf")

	#performance target:
	ls.append("analysis/dmrs/" + config['run'] + "/auc_curves.pdf")
	ls.append("analysis/dmrs/" + config['run'] + "/leave_one_out/auc_curves.pdf")

	#rms-based score using reference DMRs:
	ls.append("analysis/dmrs/" + config['run'] + "/rms/auc_curve.pdf")

	return ls

def get_multiqc_targets(wildcards):
	ls = []
	for sample in config["samples"]:
		if input_is_fq(sample):
			ls.append("analysis/fastqc/%s/%s_1_fastqc.html" % (sample, sample))
			ls.append("analysis/fastqc/%s/%s_2_fastqc.html" % (sample, sample))
			ls.append("analysis/trimgalore/%s/%s_1.fq.gz_trimming_report.txt" % (sample, sample))
			ls.append("analysis/trimgalore/%s/%s_2.fq.gz_trimming_report.txt" % (sample, sample))
			ls.append("analysis/bowtie2/%s/%s.bowtie2.log" % (sample, sample))
	return ls

def get_insert_size_file_list(wildcards):
	ls = []
	for sample in config["samples"]:
		ls.append("analysis/bowtie2/%s/%s.insert_sizes.txt" % (sample, sample))
	return ls

def get_unique_main_chr_read_file_list(wildcards):
        ls = []
        for sample in config["samples"]:
                ls.append("analysis/bowtie2/%s/%s.unique_main_chr_reads.txt" % (sample, sample))
        return ls

def get_enrichment_file_list(wildcards):
	ls = []
	for sample in config["samples"]:
		ls.append("analysis/medips/%s/%s.CpG.enrichment.txt" % (sample, sample))
	return ls

def get_coupling_set_targets():
	sample = metadata.index[0]
	if(input_is_fq(sample)):
		return "analysis/bowtie2/" + sample + "/" + sample + ".dedup.bam"
	else: # if the sample is a medips file
		return config["samples"][sample][0] 

def get_medips(wildcards):
	ls = []
	for sample in config["samples"]:
		if input_is_fq(sample):
			ls.append("analysis/medips/%s/%s.medip.rds" % (sample,sample))
		else:
			#if the input is a medips file, create a copy of it in the folder where it would be generated:
			file = config["samples"][sample][0]
			os.system("mkdir -p analysis/medips/%s" % (sample))
			os.system("cp " + file + " analysis/medips/%s/%s.medip.rds" % (sample, sample))	
	return ls

def get_DMRs_leave_one_out_targets(wildcards):                                       
        ls = []                                                                             
        for sample in list(metadata[metadata["Class"].isin(['case','control'])].index):     
                ls.append("analysis/medips/%s/%s.medip.rds" % (sample,sample))
        return ls 

def get_reference_samples(wildcards):
	ls = []
	for sample in list(metadata[metadata["Class"].isin(['reference_case','reference_control'])].index):
		ls.append("analysis/medips/%s/%s.medip.rds" % (sample,sample))
	return ls

def get_performance_leave_one_out_targets(wildcards):
	ls = []
	for sample in list(metadata[metadata["Class"].isin(['case','control'])].index):
		ls.append("analysis/dmrs/" + config['run'] + "/leave_one_out/%s/sample_prob_table.tsv" % (sample))
	return ls

# --- snakemake rules ---

# rule to run the whole workflow
rule all:
	input:
		get_targets

# copy metasheet.csv and config.yaml file to run_files directory
rule copy_run_files:
	input:
		config['metasheet']
	output:
		"run_files/metasheet.csv." + config['run']
	params:
		run = config['run']
	shell:
		"cp {input} {output} && "
		"cp config.yaml run_files/config.yaml.{params.run}"

# create links to fastq files in the attic folder
rule fastq_links:
	input:
		fq1 = get_fq1,
		fq2 = get_fq2
	output:
		fq1_ln = "attic/{sample}_1.fq.gz",
		fq2_ln = "attic/{sample}_2.fq.gz"
	shell:
		"ln -s $PWD/{input.fq1} {output.fq1_ln} && "
		"ln -s $PWD/{input.fq2} {output.fq2_ln}"

# run fastqc
rule fastqc:
	input:
		"attic/{sample}_1.fq.gz",
		"attic/{sample}_2.fq.gz"
	output:
		"analysis/fastqc/{sample}/{sample}_1_fastqc.html",
		"analysis/fastqc/{sample}/{sample}_2_fastqc.html"
	shell:
		"mkdir -p analysis/fastqc/{wildcards.sample} && "
		"fastqc -o analysis/fastqc/{wildcards.sample}/ {input}"

# perform adapter trimming
rule trim_and_fastqc:
	input:
		"attic/{sample}_1.fq.gz",
		"attic/{sample}_2.fq.gz"
	output:
		temp("analysis/trimgalore/{sample}/{sample}_1_val_1.fq.gz"),
		temp("analysis/trimgalore/{sample}/{sample}_2_val_2.fq.gz"),
		"analysis/trimgalore/{sample}/{sample}_1.fq.gz_trimming_report.txt",
		"analysis/trimgalore/{sample}/{sample}_2.fq.gz_trimming_report.txt"
	shell:
		"mkdir -p analysis/fastqc/{wildcards.sample} && "
		'trim_galore --fastqc --fastqc_args "-o analysis/fastqc/{wildcards.sample}" -o analysis/trimgalore/{wildcards.sample} --paired {input}' 

# align fastq files
### TODO should make this produce a bam file directly rather than sam file to save space/time
rule align:
	input:
		fq1 = "analysis/trimgalore/{sample}/{sample}_1_val_1.fq.gz",
		fq2 = "analysis/trimgalore/{sample}/{sample}_2_val_2.fq.gz"
	output:
		sam = temp("analysis/bowtie2/{sample}/{sample}.sam"),  ### need to test whether this really deletes sam file
		log = "analysis/bowtie2/{sample}/{sample}.bowtie2.log"
	params:
		hg19 = "ref_files/indexes/hg19"		
	shell:
		"bowtie2 -x {params.hg19} -p 4 -S {output.sam} -1 {input.fq1} -2 {input.fq2} 2> {output.log}"

# convert sam file to bam file
rule sam_to_bam:
	input:
		"analysis/bowtie2/{sample}/{sample}.sam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.bam")
	shell:
		"samtools view -S -b {input} > {output}"

# sort bam file by read names to prepare it for fixmate
rule name_sort:
	input:
		"analysis/bowtie2/{sample}/{sample}.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.namesort.bam")
	shell:
		"samtools sort -n -o {output} {input}"

# fill in mate mate coordinates and insert size fields
rule fixmate:
	input:
		"analysis/bowtie2/{sample}/{sample}.namesort.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.fixmate.bam")
	shell:
		"samtools fixmate -m {input} {output}"

# sort bam file by position to prepare for duplicate removal step
rule pos_sort:
	input:
		"analysis/bowtie2/{sample}/{sample}.fixmate.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.sorted.bam")
	shell:
		"samtools sort -o {output} {input}"

# remove duplicates reads
rule rm_dup:
	input:
		"analysis/bowtie2/{sample}/{sample}.sorted.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	shell:
		"samtools markdup -r -s {input} {output} && "
		"samtools index {output}"

# get stats on unique reads mapping to main chromosomes
rule unique_main_chr_reads:
	input:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.unique_main_chr_reads.txt"
	shell:
		""" printf "{wildcards.sample}\t%s\n" $(samtools idxstats {input} | awk '($1!="*" && $1 ~ _) {{s+=$3+$4}} END {{print s}}') > {output} """

# get insert sizes for a sampling of reads
rule insert_sizes:
	input:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.insert_sizes.txt"
	shell:
		""" samtools view -s 0.001 {input} | awk 'BEGIN{{FS="\t"}}$9>0{{print $9}}' | sort -k1,1n > {output} """

# collate qc data with multiqc
rule multiqc:
	input: 
		get_multiqc_targets
	output:
		"analysis/multiqc/multiqc.html",
		"analysis/multiqc/multiqc_data/multiqc_bowtie2.txt",
		"analysis/multiqc/multiqc_data/multiqc_fastqc.txt",
	shell:
		"multiqc analysis --ignore *val* -f -o analysis/multiqc -n multiqc.html"

# generate qc plots (alignment stats, insert sizes, read number for main chrs, CpG enrichment, duplicate reads)
rule plot_qc:
	input:
		mqc_bowtie = "analysis/multiqc/multiqc_data/multiqc_bowtie2.txt",
		mqc_fastqc = "analysis/multiqc/multiqc_data/multiqc_fastqc.txt",
		inserts = get_insert_size_file_list,
		reads = get_unique_main_chr_read_file_list,
		enrichment = get_enrichment_file_list
	output:
		alignment = "analysis/qc_plots/" + config['run'] + "/alignment_plot.pdf",
		inserts = "analysis/qc_plots/" + config['run'] + "/insert_sizes.pdf",
		reads = "analysis/qc_plots/" + config['run'] + "/main_chr_unique_reads.pdf",
		enrichment = "analysis/qc_plots/" + config['run'] + "/CpG_enrichment.pdf",
		dups = "analysis/qc_plots/" + config['run'] + "/duplication_plot.pdf"
	params:
		dir="analysis/qc_plots/" + config['run'],
		samples=",".join(metadata.index),
		meta=config['metasheet']	
	shell:
		""" mkdir -p {params.dir} &&
		Rscript scripts/plot_alignment.R "{input.mqc_bowtie}" "{output.alignment}" "{params.samples}" && 
		Rscript scripts/plot_inserts_size.R "{input.inserts}" "{output.inserts}" "{params.samples}" &&
		cat {input.reads} > {params.dir}/main_chr_unique_reads.txt &&
		Rscript scripts/plot_main_chr_unique_reads.R "{params.dir}/main_chr_unique_reads.txt" "{output.reads}" "{params.samples}" && 
		rm {params.dir}/main_chr_unique_reads.txt &&
		Rscript scripts/plot_unique_dups.R "{input.mqc_fastqc}" "{output.dups}" "{params.samples}" &&
		cat {input.enrichment} | grep GoGe > {params.dir}/CpG_enrichment.txt &&
		Rscript scripts/plot_CpG_enrichment.R "{params.dir}/CpG_enrichment.txt" "{output.enrichment}" "{params.samples}" """

# generate or load the coupling set required by MeDIPs with info about CpG content
rule coupling_set:
	input:
		get_coupling_set_targets()
	output:
		"analysis/couplingset/couplingset.rds"
	shell:
		"Rscript scripts/coupling.set.R {input} {output}"

# process bam file into medips object with read count per window
rule medips:
	input:
		bam = "analysis/bowtie2/{sample}/{sample}.dedup.bam",
		cs = "analysis/couplingset/couplingset.rds"
	output:
		"analysis/medips/{sample}/{sample}.medip.rds",
		"analysis/medips/{sample}/{sample}.CpG.enrichment.txt"
	params:
		relH = config['genome_relH'],
		GoGe = config['genome_GoGe']
	shell:
		"Rscript scripts/MeDIPS.R {wildcards.sample} {input.bam} analysis/medips/{wildcards.sample} {params.relH} {params.GoGe} && "
		"mkdir -p analysis/medips/qc && "
		"cp analysis/medips/{wildcards.sample}/*.pdf analysis/medips/qc/"

# generate pca plots, sparsitity vs depth, and other qc/analysis plots by group
rule group_plots:
	input:
		files=get_medips, #list of all medips objects, same as used for rule DMRs
		couplingset="analysis/couplingset/couplingset.rds",
		up=config['restrict_up'] if (config["restrict_up"]!="" and config["restrict_up"]!="none") else list() if (config["restrict_up"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.up",
		down=config['restrict_down'] if (config["restrict_down"]!="" and config["restrict_down"]!="none") else list() if (config["restrict_down"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.down"
	output:
		"analysis/dmrs/" + config['run'] + "/hclust.pdf",
		"analysis/dmrs/" + config['run'] + "/sparsity_depth.txt"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run']
	shell:
		""" if [ '{input.up}' = '' ]; then up='none'; else up='{input.up}'; fi
		if [ '{input.down}' = '' ]; then down='none'; else down='{input.down}'; fi
		Rscript scripts/group.medip.plots.R "{input.files}" analysis/dmrs/{params.run} {params.meta} {params.blacklist} $up $down """

# identify differentially methylated regions in a reference panel (eg, LuCaP MeDIP data)
rule reference_DMRs:
	input:
		files=get_reference_samples,
		couplingset="analysis/couplingset/couplingset.rds"
	output:
		"analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.up",
		"analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.down",
		"analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run'],
	shell:
		""" mkdir -p analysis/dmrs/{params.run}/reference
		Rscript scripts/reference_DMRs.R "{input.files}" analysis/dmrs/{params.run}/reference {params.meta} {params.blacklist} """

# identify DMRs in plasma
rule DMRs:
	input:
		files=get_medips,
		couplingset="analysis/couplingset/couplingset.rds",
		sparsitydepth="analysis/dmrs/" + config['run'] + "/sparsity_depth.txt",
		up=config['restrict_up'] if (config["restrict_up"]!="" and config["restrict_up"]!="none") else list() if (config["restrict_up"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.up",
		down=config['restrict_down'] if (config["restrict_down"]!="" and config["restrict_down"]!="none") else list() if (config["restrict_down"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.down"
	output:
		"analysis/dmrs/" + config['run'] + "/iter.{i}" + "/sample_prob_table.tsv"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run'],
	shell:
		""" if [ '{input.up}' = '' ]; then up='none'; else up='{input.up}'; fi
		if [ '{input.down}' = '' ]; then down='none'; else down='{input.down}'; fi
		Rscript scripts/DMRs.R "{input.files}" analysis/dmrs/{params.run} {params.meta} {params.blacklist} {wildcards.i} $up $down {input.sparsitydepth} """

# perform leave-one-out cross validation on plasma classification
rule DMRs_leave_one_out:
	input:
		files=get_DMRs_leave_one_out_targets,
		couplingset="analysis/couplingset/couplingset.rds",
		sparsitydepth="analysis/dmrs/" + config['run'] + "/sparsity_depth.txt",
		up=config['restrict_up'] if (config["restrict_up"]!="" and config["restrict_up"]!="none") else list() if (config["restrict_up"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.up",
		down=config['restrict_down'] if (config["restrict_down"]!="" and config["restrict_down"]!="none") else list() if (config["restrict_down"]=="none") else "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed.case.down"
	output:
		"analysis/dmrs/" + config['run'] + "/leave_one_out/{sample}/sample_prob_table.tsv"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run'],
	shell:
		""" if [ '{input.up}' = '' ]; then up='none'; else up='{input.up}'; fi
		if [ '{input.down}' = '' ]; then down='none'; else down='{input.down}'; fi
		Rscript scripts/loo_DMRs.R "{input.files}" analysis/dmrs/{params.run} {params.meta} {params.blacklist} {wildcards.sample} $up $down {input.sparsitydepth} """

# summarize and plot classification performance metrics
rule performance:
	input:
		["analysis/dmrs/" + config['run'] + "/iter.%03d/sample_prob_table.tsv" % s for s in range(1,config['iterations']+1)]
	output:
		"analysis/dmrs/" + config['run'] + "/auc_curves.pdf"
	params:
		run=config['run'],
		meta=config['metasheet'],
		iter=config['iterations']
	shell:
		"Rscript scripts/performance.plots.R analysis/dmrs/{params.run} {params.meta} {params.iter}"

# summarize and plot classification performance using leave-one-out cross-validation
rule performance_leave_one_out:
	input:
		get_performance_leave_one_out_targets
	output:
		"analysis/dmrs/" + config['run'] + "/leave_one_out/auc_curves.pdf"
	params:
		run=config['run'],
		meta=config['metasheet'],
		iter=config['iterations']
	shell:
		"Rscript scripts/performance.plots.loo.R analysis/dmrs/{params.run}/leave_one_out {params.meta} {params.iter}"

# summarize info about the DMR windows chosen for each cross-validation iteration
rule compare_DMRs:
	input:
		probs = ["analysis/dmrs/" + config['run'] + "/iter.%03d/sample_prob_table.tsv" % s for s in range(1,config['iterations']+1)],
		ref_DMRs = "analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed"
	output:
		bed = "analysis/dmrs/" + config['run'] + "/coef.summary.bed",
		hist = "analysis/dmrs/" + config['run'] + "/dmr.hist.pdf",
		plot = "analysis/dmrs/" + config['run'] + "/cf_vs_lucap.pdf"
	params:
		run=config['run']
	shell:
		""" sh scripts/plot.coefs.sh analysis/dmrs/{params.run} {output.bed} {input.ref_DMRs}
		Rscript scripts/plot.coefs.R {output.bed} analysis/dmrs/{params.run} """

# assign a score for each sample based on relative methylation score at reference (LuCaP) DMRs
rule rms:
	input:
		files=get_DMRs_leave_one_out_targets, #gets a list of case and control medips files
		ref_DMRs="analysis/dmrs/" + config['run'] + "/reference/reference_DMRs.bed"
	output:
		"analysis/dmrs/" + config['run'] + "/rms/auc_curve.pdf"
	params:
		exclude=config['rms_exclude'],
		run=config['run'],
		meta=config['metasheet']
	shell:
		""" sh scripts/subset_DMRs.sh {input.ref_DMRs} {params.exclude} analysis/dmrs/{params.run}/rms
		Rscript scripts/rms_DMRs.R "{input.files}" {params.meta} analysis/dmrs/{params.run}/rms """

