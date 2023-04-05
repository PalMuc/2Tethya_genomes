#!/bin/bash
#
#SBATCH --job-name=trimm
#SBATCH --ntasks=18
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=100G
#SBATCH --err=error_trimmomatic


java -jar ~/bin/trimmomatic.jar PE -phred33 -threads 18 \
forward.fastq.gz \
reverse.fastq.gz \
Trimmomatic_paired_forward.fastq.gz \
Trimmomatic_unpaired_forward.fastq.gz \
Trimmomatic_paired_reverse.fastq.gz \
Trimmomatic_unpaired_reverse.fastq.gz \
ILLUMINACLIP:adapter.fasta:2:30:10 SLIDINGWINDOW:10:29 HEADCROP:6 MINLEN:130 2>error_trimmomatic


fastqc Trimmomatic_paired_forward.fastq.gz
fastqc Trimmomatic_paired_reverse.fastq.gz


