## SPAdes 3.13.1

https://github.com/ablab/spades

```shell
cat spades_2.sh 
#!/bin/bash
#
#SBATCH --job-name=spades
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000

spades.py --only-assembler -1 Galaxy23-Trimmomatic_on_Taurantium_Ind01_Lib02_gDNA_Read1.fastq_R1_paired.fastq \
-2 Galaxy24-Trimmomatic_on_Taurantium_Ind01_Lib02_gDNA_Read2.fastq_R2_paired.fastq \
--nanopore ../T-citrina_nanopore/20190315_T_citrina_PCR-guppy-2.3.5.porechop.fastq --threads 224 -o assembly_spades_2 

```

```shell
perl sequence_counter.pl simplified_contigs.fasta
longest scaffold = 133249
shortest scaffold = 56
N50 = 2076
number of seuqneces = 151916
number of N = 0
N percent = 0
total length = 99966608

```