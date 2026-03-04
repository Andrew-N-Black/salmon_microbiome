#!/bin/bash
#SBATCH -A core
#SBATCH -p core
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=60G
#SBATCH --job-name=QASalmonMicro

#in working directory containing fastq files
mkdir -p fastqc_results/QC
fastqc *.fq.gz -o ./fastqc_results -t 30
cd fastqc_results/
multiqc .  -o QC/

