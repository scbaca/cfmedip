Analysis pipeline for cell-free MeDIP-seq used to classify neuroendocrine prostate cancer

Created by Sylvan Baca and Soumya Zacharia, with code adapted from Keegan Korthauer (https://github.com/kdkorthauer/cfMeDIPseq-RCC)

This pipeline is run in a conda environment specified in the file: misc/medips.env.yml

Update metasheet.csv and config.yaml with sample info as needed

To submit to cluster with slurm workload manager, use sbatch submit.sh

Prior to run, index files for bowtie need to be added to ref_files folder:

ref_files/indexes/hg19.1.bt2  ref_files/indexes/hg19.3.bt2  ref_files/indexes/hg19.rev.1.bt2
ref_files/indexes/hg19.2.bt2  ref_files/indexes/hg19.4.bt2  ref_files/indexes/hg19.rev.2.bt2

