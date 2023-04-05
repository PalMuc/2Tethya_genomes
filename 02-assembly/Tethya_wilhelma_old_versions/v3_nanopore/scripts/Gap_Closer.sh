#!/bin/bash
#
#SBATCH --job-name=GapCloser
#SBATCH --ntasks=60
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=GapCloser_wilhelma.err
#SBATCH --mem=725G


GapCloser -a scaffolded_10_x_raw.fasta \
-b Library_GapCloser.txt -o gaps_closed_scaffolded_10_x_raw.fasta -t 60
