# *Tethya minuta* assembly #

| name or version     | scaffolds | length (Mb) | N50 (kb) | longest scaffold (Mb) |
| --------            | -------- | -------- | -------- | -------- |
| v1 | 1551  | 127.3     | 180.9    | 3.190  |
| v2 | 1228  | 127.3     | 281.4    | 3.190  |
| v3 | 1043  | 139.0     | 788.3    | 4.328  |
| v4 bin 17 non-redundant <br> "Tmi-v4-no_bacteria"| 244   | 86.07     | 969.337    | 4.328             |


## metabat binning of bacteria ##
Using [metabat version 2:v2.15-25](https://bitbucket.org/berkeleylab/metabat/src/master/), and then [PRODIGAL v2.6.3](https://github.com/hyattpd/Prodigal) on all bins:

```
~/git/metabat/bin/metabat -i Tethya_minuta_V4.fasta -o TmiV4_metabat_bin -s 100000
#MetaBAT 2 (v2.15-25-g755e6aa) using minContig 2500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, maxEdges 200 and minClsSize 100000. with random seed=1681123778
~/git/metabat/bin/metabat -i Tethya_minuta_V4.fasta -o TmiV4_metabat_bin.clusters.tab -s 100000 --noBinOut --saveCls --seed 1681123778
#23 bins (135744938 bases in total) formed.
mv TmiV4_metabat_bin.17.fa TmiV4_sponge_bin.17.fa
for FILE in TmiV4_metabat_bin.*.fa ; do BASE="${FILE%.fa}" ; ~/prodigal-v2.6.3/prodigal -a $BASE.prot.fa -d $BASE.nucl.fa -f gff -o $BASE.gff -i $FILE ; done
```

bin 12 : 65% nitrite reductase to [*Paralysiella testudinis* QRQ81664.1](https://www.ncbi.nlm.nih.gov/protein/QRQ81664.1)
bin 12 : 86.64% atpB to [Gammaproteobacteria sp MXZ80730.1](https://www.ncbi.nlm.nih.gov/protein/MXZ80730.1)


| Filename            | Sequences | Total |  Mean  |  Median  | Longest  |  n50   | ID |
| --------              | --: | -------: | -------- | -------- | -------- | ------ | --- |
TmiV4_metabat_bin.17.fa | 255 | 90289220 | 354075.4 | 125264 | 4328674 | 984923 | sponge | 
| | | | | | | | |
TmiV4_metabat_bin.10.fa | 14 | 3784964 | 270354.6 | 190405 | 754858 | 487549 |  | 
TmiV4_metabat_bin.11.fa | 40 | 610267 | 15256.7 | 13338 | 44261 | 18032 |  | 
TmiV4_metabat_bin.12.fa | 8 | 4064372 | 508046.5 | 396455 | 901155 | 553638 | gamma | 
TmiV4_metabat_bin.13.fa | 54 | 10983662 | 203401.1 | 34065 | 2190704 | 1693439 |  | 
TmiV4_metabat_bin.14.fa | 30 | 1126806 | 37560.2 | 23164 | 104625 | 78867 |  | 
TmiV4_metabat_bin.15.fa | 9 | 154864 | 17207.1 | 14254 | 32313 | 16705 |  | 
TmiV4_metabat_bin.16.fa | 1 | 2071567 | 2071567.0 | 2071567 | 2071567 | 2071567 |  | 
| | | | | | | | |
TmiV4_metabat_bin.18.fa | 40 | 1075218 | 26880.5 | 21329 | 105712 | 40729 |  | 
TmiV4_metabat_bin.19.fa | 9 | 249270 | 27696.7 | 18480 | 57899 | 38181 |  | 
TmiV4_metabat_bin.1.fa | 2 | 1481924 | 740962.0 | 1467593 | 1467593 | 1467593 |  | 
TmiV4_metabat_bin.20.fa | 12 | 2938368 | 244864.0 | 137344 | 745586 | 495100 |  | 
TmiV4_metabat_bin.21.fa | 6 | 1547371 | 257895.2 | 370252 | 461187 | 417684 |  | 
TmiV4_metabat_bin.22.fa | 7 | 4423716 | 631959.4 | 743810 | 1190679 | 1190679 |  | 
TmiV4_metabat_bin.23.fa | 20 | 1933604 | 96680.2 | 22077 | 1207315 | 1207315 |  | 
TmiV4_metabat_bin.2.fa | 29 | 215338 | 7425.4 | 5994 | 17633 | 7554 |  | 
TmiV4_metabat_bin.3.fa | 3 | 705069 | 235023.0 | 140549 | 423971 | 423971 |  | 
TmiV4_metabat_bin.4.fa | 36 | 338367 | 9399.1 | 7709 | 24543 | 10938 |  | 
TmiV4_metabat_bin.5.fa | 8 | 3878909 | 484863.6 | 285074 | 1209476 | 1209476 |  | 
TmiV4_metabat_bin.6.fa | 110 | 1358553 | 12350.5 | 9125 | 59229 | 14665 |  | 
TmiV4_metabat_bin.7.fa | 8 | 335860 | 41982.5 | 43280 | 81042 | 69890 |  | 
TmiV4_metabat_bin.8.fa | 95 | 2042547 | 21500.5 | 11313 | 211448 | 30375 |  | 
TmiV4_metabat_bin.9.fa | 8 | 135102 | 16887.8 | 16482 | 29941 | 19941 |  |

`for FILE in *.gff; do Rscript ~/git/genomeGTFtools/draw_genome_annotation.R $FILE ; done`

Most bins have too many errors to correctly annotate the bacterial proteins. That is, the proteins predicted are fragmented and do not correctly blast. This was not the case for the *T. wilhelma* microbes.

## removing duplicates ##
Some scaffolds were exact duplicates, removed with [excludeAinB.py](https://bitbucket.org/wrf/sequences/src/master/excludeAinB.py) to create the non-redundant version.

```
excludeAinB.py duplicate_scaffolds_to_exclude.names TmiV4_sponge_bin.17.fa > TmiV4_sponge_bin.17.nr.fa
grep ">" TmiV4_sponge_bin.17.nr.fa | cut -c 2- > TmiV4_sponge_bin.17.nr.names
```

## wtdbg2 assembly ##

For this assembly we used all the availiable data for *Tethya minuta* we assembled and polished the reads using [wtdbg2](https://github.com/ruanjue/wtdbg2) 
and scaffolded the contigs using PE-RNA reads and [P_RNA_Scaffolder](https://github.com/CAFS-bioinformatics/P_RNA_scaffolder)

```bash
#!/bin/bash
#
#SBATCH --job-name=wtdbg2
#SBATCH --ntasks=16
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_wtdbg2_minuta
#SBATCH --mem=64G

# assembly

wtdbg2 -x ont -g 130m -i ../Tethya_minuta_nanopore_combined.fastq -t 16  -fo T_minuta_V3

# consensus sequence

wtpoa-cns -t 16 -i T_minuta_V3.ctg.lay.gz -fo T_minuta_V3.ctg.fa

# mapping

minimap2 -t 16 -x map-pb \
-a T_minuta_V3.ctg.fa ../Tethya_minuta_nanopore_combined.fastq | samtools view -Sb  >T_minuta_V3.ctg.map.bam

#polishing

samtools sort T_minuta_V3.ctg.map.bam > T_minuta_V3.ctg.map.srt.bam
samtools view T_minuta_V3.ctg.map.srt.bam | wtpoa-cns -t 16 -d T_minuta_V3.ctg.fa \
-i - -fo T_minuta_V3.ctg.2nd.fa

bwa mem -t 16  T_minuta_V3.ctg.2nd.fa \
../Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq \
../Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq \
../output_GW41566-2_S2_L001_R1_001_forward_paired.fq.gz \
../output_GW41566-2_S2_L001_R2_001_reverse_paired.fq.gz \
| samtools view -Sb - >sr.bam
samtools sort sr.bam >sr.srt.bam
samtools view sr.srt.bam | wtpoa-cns -t 16 -x sam-sr -d T_minuta_V3.ctg.2nd.fa -i - -fo T_minuta_V3.ctg.3rd.fa

bwa index T_minuta_V3.ctg.3rd.fa
bwa mem T_minuta_V3.ctg.3rd.fa  ~/T-minuta_RNA_09_2019/output_forward_paired_r1.fq.gz \
~/T-minuta_RNA_09_2019/output_reverse_paired_r4.fq.gz >RNA_T-minuta_V3.sam



``` 
## Nanopore reads mapping

* coverage 85.35 x


```bash
Longest scaffold = 3031384
Shortest scaffold = 691
N50 = 530244
Number of sequences = 1348
Number of N|n-position = 0
N|n percent = 0
Total sequence length = 124329366

```

![plot of GC content vs size](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/Tethya_minuta/figures/GC_content_vs_ln_sequence_length.png)


## P_RNA_scaffolder

```bash
#!/bin/bash
#
#SBATCH --job-name=P_RNA
#SBATCH --ntasks=32
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_P_RNA_minuta
#SBATCH --mem=200G



sh ~/program/P_RNA_scaffolder/P_RNA_scaffolder.sh \
-d ~/program/P_RNA_scaffolder  -i RNA_T-minuta_V3.sam -j T_minuta_V3.ctg.3rd.fa \
-o P_RNA_results_1000 -t 32  -F ~/T-minuta_RNA_09_2019/output_forward_paired_r1.fq.gz \
-R ~/T-minuta_RNA_09_2019/output_reverse_paired_r4.fq.gz -b no  -e 100000 -s yes 2>P_RNA_err

```
```bash
Longest scaffold = 4328638
Shortest scaffold = 691
N50 = 654661
Number of sequences = 1250
Number of N|n-position = 9800
N|n percent = 0.00788166779243155
Total sequence length = 124339166


```
## Short read mapping 

```bash
42728251 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
39291 + 0 supplementary
0 + 0 duplicates
36836636 + 0 mapped (86.21% : N/A)
42688960 + 0 paired in sequencing
21344480 + 0 read1
21344480 + 0 read2
35421352 + 0 properly paired (82.98% : N/A)
35724766 + 0 with itself and mate mapped
1072579 + 0 singletons (2.51% : N/A)
207308 + 0 with mate mapped to a different chr
129320 + 0 with mate mapped to a different chr (mapQ>=5)

```
## Long read mapping 

```bash
SN      raw total sequences:    6937646
SN      filtered sequences:     0
SN      sequences:      6937646
SN      is sorted:      1
SN      1st fragments:  6937646
SN      last fragments: 0
SN      reads mapped:   5160415
SN      reads mapped and paired:        0       # paired-end technology bit set + both mates mapped
SN      reads unmapped: 1777231
SN      reads properly paired:  0       # proper-pair bit set
SN      reads paired:   0       # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      10226   # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 1286615
SN      total length:   14611841287     # ignores clipping
SN      bases mapped:   11985348431     # ignores clipping
SN      bases mapped (cigar):   11055883132     # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     1351213865      # from NM fields
SN      error rate:     1.222167e-01    # mismatches / bases mapped (cigar)
SN      average length: 2106
SN      maximum length: 2581795
SN      average quality:        16.4
SN      insert size average:    0.0
SN      insert size standard deviation: 0.0
SN      inward oriented pairs:  0
SN      outward oriented pairs: 0
SN      pairs with other orientation:   0
SN      pairs on different chromosomes: 0
```

## SSPACE_Standard short read scaffolding ##
```bash
#!/bin/bash
#
#SBATCH --job-name=sspace
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --time=48:00:00
#SBATCH --mem=55000

#scaffolding with short reads

perl SSPACE_Standard_v3.0.pl \
-l Library_minuta_PE_2.txt \
-s Tethya_minuta_P_RNA.fasta \
-x 1 \
-o 5 \
-b output_minuta_3 \
-T 28 \
>logfile_minuta_3_sspace 2>logfile_3_minuta_sspace
```
```bash
                Longest scaffold = 4328674
                Shortest scaffold = 1040
                N50 = 654661
                Number of sequences = 1238
                Number of N|n-position = 9867
                N|n percent = 0.00793276941341508
                Total sequence length = 124382791

```

## LRScaf ##
```bash
#!/bin/bash
#
#SBATCH --job-name=LRScaf
#SBATCH --ntasks=32
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_LRScaf_minuta
#SBATCH --mem=200G


minimap2  -x map-ont -o mapped_longreads  output_minuta_3.final.scaffolds.fasta   Tethya_minuta_nanopore_combined.fastq

java -jar ~/program/lrscaf/target/LRScaf-1.1.9.jar  -c output_minuta_3.final.scaffolds.fasta -a mapped_longreads \
-t mm -o scaffolded_minuta



```
```bash
9334911 + 0 in total (QC-passed reads + QC-failed reads)
1461609 + 0 secondary
935656 + 0 supplementary
0 + 0 duplicates
7737278 + 0 mapped (82.89% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
```bash
Longest scaffold = 4328674
Shortest scaffold = 2616
N50 = 788386
Number of sequences = 1043
Number of N|n-position = 887046
N|n percent = 0.637792656422024
Total sequence length = 139080623

```
![cumulative length plot](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/Tethya_minuta/figures/cumulative_length_plot.png)
![histogram of contig GC percent](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/Tethya_minuta/figures/GC_histogram.png)

## BUSCO ##

```bash
#!/bin/bash
#
#SBATCH --job-name=busco
#SBATCH --ntasks=4
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=err_busco
#SBATCH --mem=10G

export AUGUSTUS_CONFIG_PATH=~/program/augustus-3.3.1/config

python ~/program/busco-master/scripts/run_BUSCO.py --in ../scaffolds.fasta  -o BUSCO -l ~/program/busco-master/metazoa_odb9 \
-c 4 -f -m geno >logfile 2>error

```
```bash
INFO    Results:
INFO    C:71.2%[S:62.9%,D:8.3%],F:10.6%,M:18.2%,n:978
INFO    696 Complete BUSCOs (C)
INFO    615 Complete and single-copy BUSCOs (S)
INFO    81 Complete and duplicated BUSCOs (D)
INFO    104 Fragmented BUSCOs (F)
INFO    178 Missing BUSCOs (M)
INFO    978 Total BUSCO groups searched
INFO    BUSCO analysis done. Total running time: 4353.883141994476 seconds
INFO    Results written in /home/fdeiste/Tethya_minuta_all/PromethION/LRScaf_2/scaffolded_minuta/BUSCO/run_BUSCO/
```
## Hisat2 RNA ##
```bash
#!/bin/bash
#
#SBATCH --job-name=Hisat2
#SBATCH --ntasks=16
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=Hisat2.err
#SBATCH --mem=200G


# build index for query sequence for Hisat2

#hisat2-build ../scaffolds.fasta \
#scaffolds.fasta

# run HISAT v. 2.1.0

hisat2 -q --no-mixed --no-unal --no-spliced-alignment --threads 16 -x scaffolds.fasta \
-1 ../../../../output_forward_paired_r1.fq.gz  -2 ../../../../output_reverse_paired_r4.fq.gz \
-S scaffolds.sam 2>errorlog.txt

samtools stats -@16  scaffolds.sam >samstats.txt
samtools flagstat -@16  scaffolds.sam >flagstats.txt
```
## samtools stats ##
```bash
SN      raw total sequences:    165241566
SN      filtered sequences:     0
SN      sequences:      165241566
SN      is sorted:      0
SN      1st fragments:  82620783
SN      last fragments: 82620783
SN      reads mapped:   165241566
SN      reads mapped and paired:        165241566       # paired-end technology bit set + both mates mapped
SN      reads unmapped: 0
SN      reads properly paired:  164417228       # proper-pair bit set
SN      reads paired:   165241566       # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      4564180 # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 34954762
SN      total length:   7270628904      # ignores clipping
SN      bases mapped:   7270628904      # ignores clipping
SN      bases mapped (cigar):   7002843460      # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     29933494        # from NM fields
SN      error rate:     4.274477e-03    # mismatches / bases mapped (cigar)
SN      average length: 44
SN      maximum length: 44
SN      average quality:        36.8
SN      insert size average:    338.4
SN      insert size standard deviation: 579.0
SN      inward oriented pairs:  82115327
SN      outward oriented pairs: 339938
SN      pairs with other orientation:   9025
SN      pairs on different chromosomes: 63049

```
## BRAKER2 ##
```bash
#!/bin/bash
#
#SBATCH --job-name=braker
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --time=48:00:00
#SBATCH --mem=90000


export PATH=$PATH:/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/hisat2-2.1.0

# build index for query sequence for Hisat2

hisat2-build Tethya_minuta_all_LRScaf_2.fasta \
scaffolds.fasta

# run HISAT v. 2.1.0

hisat2 -q --no-mixed --no-unal --threads 28 -x scaffolds.fasta \
-1 output_forward_paired_r1.fq.gz  -2 output_reverse_paired_r4.fq.gz \
-S scaffolds.sam 2>errorlog.txt

samtools stats -@28  scaffolds.sam >samstats.txt
samtools flagstat -@28  scaffolds.sam >flagstats.txt


samtools sort -n scaffolds.bam >sort_scaffolds.bam 
 
braker.pl --GENEMARK_PATH=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/gmes_linux_64  \
--genome=renamed_Tethya_minuta_all_LRScaf_2.fasta --species=Tethya_minuta --bam=sort_scaffolds.bam  

```
## AUGUSTUS ##
```bash
```
## OrthoFinder ##
```bash
```

## Illumina DNA reads ##

* 100 bp :  42,000,000 reads
* 150 bp :   8,000,000 reads

## Illumina RNA reads ##

* 50 bp  : 274,000,000 reads

## Nanopore 

* minION : 2.17 GBases
 
* PromethION : 12.73 GBases

## Scaffolds statistic *Tethya minuta* V3

![plot of T minuta GC versus coverage](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/Tethya_minuta/figures/Tethya_minuta_all_scaffolded_LRS_P_RNAscaffolds.fasta.png)


