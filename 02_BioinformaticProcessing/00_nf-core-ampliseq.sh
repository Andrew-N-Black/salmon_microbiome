#!/bin/bash
#SBATCH -A core
#SBATCH -p core
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=100G
#SBATCH --job-name=SM_nfcore

export SINGULARITY_TMPDIR=/scratch
nextflow run /local/cqls/software/nextflow/assets/nf-core-ampliseq_2.14.0/2_14_0/ -profile singularity -c nextflow.config -params-file  nf-params.json -resume
