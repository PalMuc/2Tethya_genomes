#!/bin/bash
#
#SBATCH --job-name=trimm
#SBATCH --ntasks=1
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --err=error_trimmomatic


gzip forward.fastq 
gzip reverse.fastq


