## P_RNA_scaffolder

We used P_RNA_scaffolder to scaffold *Tethya minuta* V1. before scaffolding, the statistics were as follows

```shell
                Longest scaffold = 3190713
                Shortest scaffold = 3434
                N50 = 182355
                Number of sequences = 1550
                Number of N|n-position = 1277938
                N|n percent = 1.00377971244715
                Total sequence length = 127312595

```
```shell

#!/bin/bash
#SBATCH --job-name=P_RNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000

sh /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/P_RNA_scaffolder/P_RNA_scaffolder.sh \
-d /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/P_RNA_scaffolder -i RNA_T-minuta.sam -j \
standard_output.final.scaffolds.fasta \
-o P_RNA_results_100000 -t 28 -e 100000 -F link_Trimmomatic_paired_GW41610_TGCACGAT_1.fastq \
-R link_Trimmomatic_paired_GW41610_TGCACGAT_2.fastq -b no  -s yes 2>P_RNA_err

```

```shell
perl sequence_counter.pl T-minuta_P_RNA_100000_scaffold.fasta

                Longest scaffold = 3190713
                Shortest scaffold = 3434
                N50 = 281420
                Number of sequences = 1228
                Number of N|n-position = 1310138
                N|n percent = 1.02881158197318
                Total sequence length = 127344795
```

