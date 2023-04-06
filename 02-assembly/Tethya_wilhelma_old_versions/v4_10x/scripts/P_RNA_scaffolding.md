## RNA scaffolding

We used P_RNA_scaffolder on the LRZ cluster  (https://github.com/CAFS-bioinformatics/P_RNA_scaffolder/blob/master/P_RNA_scaffolder.sh) to scaffold *Tethya wilhelma* V3 to become *Tethya wilhelma* V4

To use python version 2 for P_RNA scaffolder, we crated an environemt using python 2
```bash
conda create --name py2 python=2.7
conda activate py2
```
RNA reads were mapped to the genome draft using BWA (https://github.com/lh3/bwa.git)

```bash
bwa index Output_gaps_closed_scaffolded_10_x_raw_without_species_name.fasta
bwa mem Output_gaps_closed_scaffolded_10_x_raw_without_species_name.fasta Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq \
Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq >RNA_T-wilhelma.sam
```
P_RNA_scaffolder
```bash
#!/bin/bash
#
#SBATCH --job-name=P_RNA
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000


sh /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/P_RNA_scaffolder/P_RNA_scaffolder.sh \
-d /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/P_RNA_scaffolder -i RNA_T-wilhelma.sam -j \
Output_gaps_closed_scaffolded_10_x_raw_without_species_name.fasta \
-o P_RNA_results_1000 -t 28 -e 10 -F Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq \
-R Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq -b no -e 1000 -s yes 2>P_RNA_err
```

```bash
#!/bin/bash
#
#SBATCH --job-name=bwa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000

# mapping using bwa

bwa mem -M P_RNA_scaffold.fasta  ../Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq \
../Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq |samtools view -Sb - > P_RNA_1000_RNA_T-wilhelma.bam

samtools flagstat P_RNA_1000_RNA_T-wilhelma.bam >flagstat.txt
```

```bash 
cat flagstat.txt 
237305498 + 0 in total (QC-passed reads + QC-failed reads)
39365382 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
224608585 + 0 mapped (94.65% : N/A)
197940116 + 0 paired in sequencing
98970058 + 0 read1
98970058 + 0 read2
169714632 + 0 properly paired (85.74% : N/A)
183839570 + 0 with itself and mate mapped
1403633 + 0 singletons (0.71% : N/A)
1921096 + 0 with mate mapped to a different chr
1011724 + 0 with mate mapped to a different chr (mapQ>=5)
```

```shell
# BUSCO version is: 3.1.0 
# The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# To reproduce this run: python ../scripts/run_BUSCO.py -i Tethya_wilhelma_V4.fasta -o OUTPUT -l /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/busco-master/metazoa_odb9/ -m genome -c 28 -sp fly
#
# Summarized benchmarking in BUSCO notation for file Tethya_wilhelma_V4.fasta
# BUSCO was run in mode: genome

        C:82.7%[S:79.0%,D:3.7%],F:3.4%,M:13.9%,n:978

        809     Complete BUSCOs (C)
        773     Complete and single-copy BUSCOs (S)
        36      Complete and duplicated BUSCOs (D)
        33      Fragmented BUSCOs (F)
        136     Missing BUSCOs (M)
        978     Total BUSCO groups searched
```