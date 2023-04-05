#!/bin/bash
#
#SBATCH --job-name=10x
#SBATCH --ntasks=32
#SBATCH --mail-user f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --mem=125G
#SBATCH --err=error_mkoutput_pseudohap2.log



supernova mkoutput --outprefix=Twilhelma_10x_pseudohap2 --asmdir=Twilhelma_final/outs/assembly --style=pseudohap2 --nozip --headers=full
