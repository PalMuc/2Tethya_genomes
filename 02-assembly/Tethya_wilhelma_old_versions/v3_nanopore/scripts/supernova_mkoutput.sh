#!/bin/bash
#
#SBATCH --job-name=10x
#SBATCH --ntasks=32
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=125G
#SBATCH --err=error_mkoutput.log



supernova mkoutput --outprefix=Twilhelma_10x_raw --asmdir=Twilhelma_final/outs/assembly --style=raw --nozip --headers=full
