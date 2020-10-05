configfile: "config.yaml"


import pandas as pd
import yaml
import os

metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#')

def input_is_fq(sample):
	first_file = config["samples"][sample][0]
	if first_file.endswith('.fq.gz'):
		return True
	else:
		return False

# Get classes labels for each sample
def get_Class(cls):
    comp = metadata["Class"]
    return metadata[comp == cls].index

def get_fq1(wildcards):
	return config["samples"][wildcards.sample][0]

def get_fq2(wildcards):
	if(len(config["samples"][wildcards.sample])) == 1:
		return config["samples"][wildcards.sample][0]
		print("only one file provided for this sample, assuming it is a medips object")
	else:
		return config["samples"][wildcards.sample][1]

# create targets for "all rule"
def get_targets(wildcards):
	ls = []
	
	#multiqc target:
	ls.append("analysis/multiqc/multiqc.html")

	#plot_qc target:
	ls.append("analysis/qc_plots/" + config['run'] + "/main_chr_unique_reads.pdf")
	ls.append("analysis/qc_plots/" + config['run'] + "/CpG_enrichment.pdf")

	#DMRs targets:
#	ls.append("analysis/dmrs/" + config['run'] + "/case.control.diff.rds") #eventually will make this accept labels other than case/control

	#group_plots targets:
	ls.append("analysis/dmrs/" + config['run'] + "/hclust.pdf")

	#temp - can remove this once a step plotting the data is added
	for sample in config["samples"]:
        	ls.append("analysis/bowtie2/%s/%s.unique_main_chr_reads.txt" % (sample, sample))
        	ls.append("analysis/bowtie2/%s/%s.insert_sizes.txt" % (sample, sample))

	#performance target:
	ls.append("analysis/dmrs/" + config['run'] + "/auc_curves.pdf")
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

rule all:
	input:
		get_targets

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

rule fastq_links:
	input:
		fq1 = get_fq1,
		fq2 = get_fq2
	output:
		fq1_ln = "attic/{sample}_1.fq.gz",
		fq2_ln = "attic/{sample}_2.fq.gz"
	shell:
		"ln -s {input.fq1} {output.fq1_ln} && "
		"ln -s {input.fq2} {output.fq2_ln}"

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
	
rule sam_to_bam:
	input:
		"analysis/bowtie2/{sample}/{sample}.sam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.bam")
	shell:
		"samtools view -S -b {input} > {output}"

rule name_sort:
	input:
		"analysis/bowtie2/{sample}/{sample}.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.namesort.bam")
	shell:
		"samtools sort -n -o {output} {input}"

rule fixmate:
	input:
		"analysis/bowtie2/{sample}/{sample}.namesort.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.fixmate.bam")
	shell:
		"samtools fixmate -m {input} {output}"

rule pos_sort:
	input:
		"analysis/bowtie2/{sample}/{sample}.fixmate.bam"
	output:
		temp("analysis/bowtie2/{sample}/{sample}.sorted.bam")
	shell:
		"samtools sort -o {output} {input}"

rule rm_dup:
	input:
		"analysis/bowtie2/{sample}/{sample}.sorted.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	shell:
		"samtools markdup -r -s {input} {output} && "
		"samtools index {output}"

rule unique_main_chr_reads:
	input:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.unique_main_chr_reads.txt"
	shell:
		""" printf "{wildcards.sample}\t%s\n" $(samtools idxstats {input} | awk '($1!="*" && $1 ~ _) {{s+=$3+$4}} END {{print s}}') > {output} """


rule insert_sizes:
	input:
		"analysis/bowtie2/{sample}/{sample}.dedup.bam"
	output:
		"analysis/bowtie2/{sample}/{sample}.insert_sizes.txt"
	shell:
		""" samtools view -s 0.001 {input} | awk 'BEGIN{{FS="\t"}}$9>0{{print $9}}' | sort -k1,1n > {output} """


rule multiqc:
	input: 
		get_multiqc_targets
	output:
		"analysis/multiqc/multiqc.html",
		"analysis/multiqc/multiqc_data/multiqc_general_stats.txt"
	shell:
		"multiqc analysis --ignore *val* -f -o analysis/multiqc -n multiqc.html"

#TODO: update these to plot the qc metrics by batch and/or by case vs control group
rule plot_qc:
	input:
		mqc = "analysis/multiqc/multiqc_data/multiqc_general_stats.txt",
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
		Rscript scripts/plot_alignment.R "{input.mqc}" "{output.alignment}" "{params.samples}" && 
		Rscript scripts/plot_inserts_size.R "{input.inserts}" "{output.inserts}" "{params.samples}" &&
		cat {input.reads} > {params.dir}/main_chr_unique_reads.txt &&
		Rscript scripts/plot_main_chr_unique_reads.R "{params.dir}/main_chr_unique_reads.txt" "{output.reads}" "{params.samples}" && 
		rm {params.dir}/main_chr_unique_reads.txt &&
		Rscript scripts/plot_unique_dups.R "{input.mqc}" "{output.dups}" "{params.samples}" &&
		cat {input.enrichment} | grep GoGe > {params.dir}/CpG_enrichment.txt &&
		Rscript scripts/plot_CpG_enrichment.R "{params.dir}/CpG_enrichment.txt" "{output.enrichment}" "{params.samples}" """

rule coupling_set:
	input:
		get_coupling_set_targets()
	output:
		"analysis/couplingset/couplingset.rds"
	shell:
		"Rscript scripts/coupling.set.R {input} {output}"

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

rule group_plots:
	input:
		files=get_medips, #list of all medips objects, same as used for rule DMRs
		couplingset="analysis/couplingset/couplingset.rds"
	output:
		"analysis/dmrs/" + config['run'] + "/hclust.pdf"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run'],
                up=config['restrict_up'],                                                                                                 
                down=config['restrict_down']
	shell:
		""" Rscript scripts/group.medip.plots.R "{input.files}" analysis/dmrs/{params.run} {params.meta} {params.blacklist} {params.up} {params.down} """

rule DMRs:
	input:
		files=get_medips,
		couplingset="analysis/couplingset/couplingset.rds"
	output:
		"analysis/dmrs/" + config['run'] + "/iter.{i}" + "/performance.tsv"
	params:
		meta=config['metasheet'],
		blacklist=config['blacklist'],
		run=config['run'],
		up=config['restrict_up'],
		down=config['restrict_down'],
	shell:
		""" Rscript scripts/DMRs.R "{input.files}" analysis/dmrs/{params.run} {params.meta} {params.blacklist} {wildcards.i} {params.up} {params.down} """

rule performance:
	input:
		["analysis/dmrs/" + config['run'] + "/iter.%03d/performance.tsv" % s for s in range(1,config['iterations']+1)]
	output:
		"analysis/dmrs/" + config['run'] + "/auc_curves.pdf"
	params:
		run=config['run'],
		meta=config['metasheet'],
		iter=config['iterations']
	shell:
		"Rscript scripts/performance.plots.R analysis/dmrs/{params.run} {params.meta} {params.iter}"
