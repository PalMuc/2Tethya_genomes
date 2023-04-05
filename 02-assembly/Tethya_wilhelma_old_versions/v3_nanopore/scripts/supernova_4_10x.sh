#!/bin/bash
#
#SBATCH --job-name=10x
#SBATCH --ntasks=80
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=755G
#SBATCH --err=error2.log



supernova run --id=Twilhelma_final --fastqs=~/Tethya_wilhelma_10x/all_data --sample=TW --maxreads="all" --localcores=80 --localmem=755 --accept-extreme-coverage

