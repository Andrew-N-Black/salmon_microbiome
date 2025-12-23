#!/bin/bash
#SBATCH -A core
#SBATCH -p core
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=100G
#SBATCH --job-name=SM_decon

#bwa index ../salmon_ref/GCF_018296145.1_Otsh_v2.0_genomic.fna
for i in $(ls -1 ../../ALab_0006/*.fastq.gz | sed 's/.\{14\}$//' | uniq |  sed 's|../../ALab_0006/||g')
do

echo "Aligning reads and removing salmon dna from sample $i"

bwa mem -t 30 ../salmon_ref/GCF_018296145.1_Otsh_v2.0_genomic.fna ../../ALab_0006/${i}1_001.fastq.gz ../../ALab_0006/${i}2_001.fastq.gz > mapping/${i}.sam

samtools view -Sb -f 12 -F 256 mapping/${i}.sam | bedtools bamtofastq -i stdin -fq cleaned/${i}1_001.fastq -fq2 cleaned/${i}2_001.fastq

samtools flagstat  mapping/${i}.sam -@ 30 -O mapping/${i}.txt

done
