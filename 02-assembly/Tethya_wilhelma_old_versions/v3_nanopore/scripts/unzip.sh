#!/bin/bash
#
#SBATCH --job-name=gunzip
#SBATCH --ntasks=32
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=125G
#SBATCH --err=error_gunzip.log

gunzip -c TW_S0_L001_R1_001.fastq.gz >TW_S0_L001_R1_001.fastq
gunzip -c TW_S0_L002_R1_001.fastq.gz >TW_S0_L002_R1_001.fastq
gunzip -c TW_S0_L003_R1_001.fastq.gz >TW_S0_L003_R1_001.fastq
gunzip -c TW_S0_L004_R1_001.fastq.gz >TW_S0_L004_R1_001.fastq

cat TW_S0_L001_R1_001.fastq TW_S0_L002_R1_001.fastq TW_S0_L003_R1_001.fastq TW_S0_L004_R1_001.fastq >forward.fastq 

