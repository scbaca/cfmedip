#!/bin/bash
#
#SBATCH --job-name=medips
#SBATCH --output=out.medips.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium
#SBATCH --mem=8G
#SBATCH -t 2-00:00

srun snakemake -s medips.Snakefile -j 100 -pr --rerun-incomplete --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads}"
