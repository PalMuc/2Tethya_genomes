```perl
#!/usr/bin/perl -w 

open(OUT,">>Tethya_citrina_nanopore_combined.fastq");

open(DIR,"<20190315_T_citrina_PCR-guppy-2.3.5.porechop.fastq");

while(<DIR>)
        {
                my$line=$_;
                chomp$line;
                print OUT "$line\n";
        }
close DIR;
print OUT
open(DIR,"<T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.fastq");

while(<DIR>)
        {
                my$line=$_;
                chomp $line;
                print OUT "$line\n";
        }
```

```bash
#!/bin/bash
#
#SBATCH --job-name=SPAdes
#SBATCH --ntasks=16
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_spades
#SBATCH --mem=20G


~/program/SPAdes-3.13.0-Linux/bin/spades.py --pe1-1 C840RACXX_Taurantium_15s015355-1-1_Vargas-R_lane515s015355_1_sequence.fq.gz \
--pe1-2 C840RACXX_Taurantium_15s015355-1-1_Vargas-R_lane515s015355_2_sequence.fq.gz --pe2-1 C840RACXX_Taurantium_15s015356-1-1_Vargas-R_lane515s015356_1_sequence.fq.gz \
--pe2-2 C840RACXX_Taurantium_15s015356-1-1_Vargas-R_lane515s015356_2_sequence.fq.gz --nanopore Tethya_citrina_nanopore_combined.fastq \
-t 16 -m 20 -o spades_output

```