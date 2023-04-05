#!/bin/bash
#
#SBATCH --job-name=gunzip
#SBATCH --ntasks=8
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=77G
#SBATCH --err=error_gunzip_reverse.log

gunzip -c TW_S0_L001_R2_001.fastq.gz >TW_S0_L001_R2_001.fastq
gunzip -c TW_S0_L002_R2_001.fastq.gz >TW_S0_L002_R2_001.fastq
gunzip -c TW_S0_L003_R2_001.fastq.gz >TW_S0_L003_R2_001.fastq
gunzip -c TW_S0_L004_R2_001.fastq.gz >TW_S0_L004_R2_001.fastq

cat TW_S0_L001_R2_001.fastq TW_S0_L002_R2_001.fastq TW_S0_L003_R2_001.fastq TW_S0_L004_R2_001.fastq >reverse.fastq 

