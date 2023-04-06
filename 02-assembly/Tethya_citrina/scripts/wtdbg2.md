```bash
#!/bin/bash
#
#SBATCH --job-name=wtdbg2
#SBATCH --ntasks=16
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_spades
#SBATCH --mem=64G


wtdbg2 -x ont -g 100000000 -i ../Tethya_citrina_nanopore_combined.fastq -t 16  -fo T_citrina_PCR-guppy-2.3.5

wtpoa-cns -t 16 -i T_citrina_PCR-guppy-2.3.5.ctg.lay.gz -fo T_citrina_PCR-guppy-2.3.5.ctg.fa

minimap2 -t 16 -x map-pb \
-a T_citrina_PCR-guppy-2.3.5.ctg.fa ../Tethya_citrina_nanopore_combined.fastq | samtools view -Sb  >T_citrina_PCR-guppy-2.3.5.ctg.map.bam

samtools sort T_citrina_PCR-guppy-2.3.5.ctg.map.bam > T_citrina_PCR-guppy-2.3.5.ctg.map.srt.bam
samtools view T_citrina_PCR-guppy-2.3.5.ctg.map.srt.bam | wtpoa-cns -t 16 -d T_citrina_PCR-guppy-2.3.5.ctg.fa \
-i - -fo T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa

bwa mem -t 16  T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa \
C840RACXX_Taurantium_15s015355-1-1_Vargas-R_lane515s015355_1_sequence.fq.gz \
C840RACXX_Taurantium_15s015355-1-1_Vargas-R_lane515s015355_2_sequence.fq.gz \
C840RACXX_Taurantium_15s015356-1-1_Vargas-R_lane515s015356_1_sequence.fq.gz \
C840RACXX_Taurantium_15s015356-1-1_Vargas-R_lane515s015356_2_sequence.fq.gz \
| samtools view -Sb - >sr.bam
samtools sort sr.bam >sr.srt.bam
samtools view sr.srt.bam | wtpoa-cns -t 16 -x sam-sr -d T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa -i - -fo T_citrina_PCR-guppy-2.3.5.ctg.3rd.fa


```
