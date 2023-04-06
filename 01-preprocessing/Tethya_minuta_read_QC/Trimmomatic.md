```bash

#!/bin/bash
#
#SBATCH --job-name=hm.batchA1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=1:00:00
#SBATCH --mem=55000

trimmomatic PE -phred33 -threads 56 \
L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq \
L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq \
Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq \
Trimmomatic_unpaired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq \
Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq \
Trimmomatic_unpaired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq \
-trimlog trimmomatic_log -summary summary_trimmomatic \
ILLUMINACLIP:adapter.fasta:2:30:10 MINLEN:36 \
SLIDINGWINDOW:10:29 HEADCROP:7 

fastqc Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq

fastqc Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq
```