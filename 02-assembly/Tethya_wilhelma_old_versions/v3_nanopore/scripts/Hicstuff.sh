#!/bin/bash
#
#SBATCH --job-name=Hi-C 
#SBATCH --ntasks=32
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_Hi-C
#SBATCH --mem=125G

export PYTHONPATH=/home/fdeiste/.local/lib64/python3.6/

hicstuff pipeline -t 32 -a minimap2   --plot  -e DpnII -o out/ --iterative --no-cleanup -g gaps_closed_scaffolded_10_x_raw.fasta  output_DTG-HiC-396_R1_001_forward_paired.fq.gz output_DTG-HiC-396_R2_001_reverse_paired.fq.gz


