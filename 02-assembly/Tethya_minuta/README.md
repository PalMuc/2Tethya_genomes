# *Tethya minuta* assembly #

| name or version     | scaffolds | length (Mb) | N50 (kb) | longest scaffold (Mb) |
| --------            | -------- | -------- | -------- | -------- |
| "final" V1          | 1712  |  93.972     | 99.8    | 1.766                 |
| scaffolded nanopore | 1551  | 127.312     | 180.9    | 3.190                 |
| gap closed nanopore | 1551  | 127.312     | 180.9    | 3.190                 |
| V4                  | 1043  | 139.080     | 788.3    | 4.328                 |

## wtdbg2 assembly ##

For this assembly we used all the availiable data for *Tethya minuta* we assembled and polished the reads using wtdbg2 (https://github.com/ruanjue/wtdbg2) 
and scaffolded the contigs using PE-RNA reads and P_RNA_Scaffolder (https://github.com/CAFS-bioinformatics/P_RNA_scaffolder)

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
## SSPACE_Standard short read scaffolding
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

## LRScaf
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

## BUSCO

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
## Hisat2 RNA
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
## samtools stats
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
## BRAKER2
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
## AUGUSTUS
```bash
```
## OrthoFinder
```bash
```

## Illumina DNA reads

* 100 bp :  42,000,000 reads
* 150 bp :   8,000,000 reads

## Illumina RNA reads

* 50 bp  : 274,000,000 reads

## Nanopore 

* minION : 2.17 GBases
 
* PromethION : 12.73 GBases

## Scaffolds statistic *Tethya minuta* V3

![plot of T minuta GC versus coverage](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/Tethya_minuta/figures/Tethya_minuta_all_scaffolded_LRS_P_RNAscaffolds.fasta.png)

| name | length |GC content |
| -------- | -------- | -------- |
| Scaffolds_1040 | 4328674 | 39.1193519681307 |
| Scaffolds_23 | 3614561 | 38.7681455615635 |
| Scaffolds_22 | 2416165 | 38.9230609662497 |
| Scaffolds_940 | 2190704 | 69.4292616170312 |
| Scaffolds_106 | 2148092 | 39.0830046342817 |
| Scaffolds_101 | 2071567 | 53.2213618105979 |
| Scaffolds_146 | 2027610 | 38.7311240708253 |
| Scaffolds_181 | 1992245 | 38.718926063727 |
| Scaffolds_150 | 1967258 | 38.5910072576692 |
| Scaffolds_675 | 1926824 | 39.1152713163811 |
| Scaffolds_644 | 1777800 | 39.3404889858376 |
| Scaffolds_184 | 1696831 | 69.2516314491145 |
| Scaffolds_14 | 1693439 | 70.0344803839386 |
| Scaffolds_375 | 1631807 | 39.3471425305994 |
| Scaffolds_710 | 1587595 | 38.6417963752766 |
| Scaffolds_249 | 1577474 | 38.7617450132427 |
| Scaffolds_647 | 1577474 | 38.7617307648421 |
| Scaffolds_123 | 1519745 | 38.7933789140384 |
| Scaffolds_746 | 1482276 | 38.9764416979248 |
| Scaffolds_1018 | 1467593 | 38.3682305155979 |
| Scaffolds_520 | 1367625 | 38.7306626733751 |
| Scaffolds_8 | 1273932 | 38.1255416284649 |
| Scaffolds_744 | 1224589 | 39.6561790834654 |
| Scaffolds_136 | 1209476 | 62.2480735522704 |
| Scaffolds_348 | 1209476 | 62.2480735522704 |
| Scaffolds_822 | 1207315 | 45.1829138494695 |
| Scaffolds_96 | 1197368 | 39.0546964518963 |
| Scaffolds_393 | 1197368 | 39.0546781696552 |
| Scaffolds_236 | 1190679 | 59.4809869629448 |
| Scaffolds_676 | 1190679 | 59.4810028882631 |
| Scaffolds_200 | 1150242 | 38.2033008088922 |
| Scaffolds_135 | 1077559 | 38.4520885673823 |
| Scaffolds_28 | 1044500 | 38.3170125946712 |
| Scaffolds_58 | 1043729 | 39.0299153455678 |
| Scaffolds_174 | 1043281 | 38.3282388383593 |
| Scaffolds_210 | 1030245 | 38.0146751136154 |
| Scaffolds_224 | 984923 | 38.2571561048468 |
| Scaffolds_137 | 969337 | 38.7898346273877 |
| Scaffolds_44 | 938770 | 37.8690717893764 |
| Scaffolds_396 | 916768 | 38.20794746774 |
| Scaffolds_81 | 901155 | 63.8253933554308 |
| Scaffolds_165 | 901155 | 63.8253933554308 |
| Scaffolds_211 | 881861 | 39.2965020085771 |
| Scaffolds_193 | 876994 | 38.3284307222164 |
| Scaffolds_258 | 861833 | 38.1609904238856 |
| Scaffolds_195 | 858504 | 38.0205345909378 |
| Scaffolds_814 | 851675 | 38.3452881784996 |
| Scaffolds_53 | 788386 | 38.6895792427079 |
| Scaffolds_251 | 773658 | 39.5190664657176 |
| Scaffolds_609 | 773658 | 39.5190393713001 |
| Scaffolds_10 | 773467 | 37.5701623186676 |
| Scaffolds_877 | 771878 | 58.7363320723429 |
| Scaffolds_253 | 763023 | 39.153869028654 |
| Scaffolds_69 | 754858 | 58.7152123024904 |
| Scaffolds_103 | 745586 | 59.659490227847 |
| Scaffolds_0 | 743810 | 59.8127335770855 |
| Scaffolds_516 | 735519 | 37.9347428557444 |
| Scaffolds_259 | 734429 | 69.9681250706327 |
| Scaffolds_961 | 734429 | 69.968179447763 |
| Scaffolds_40 | 701041 | 38.8849196443174 |
| Scaffolds_177 | 696990 | 38.0516242543252 |
| Scaffolds_160 | 675768 | 38.2084238904526 |
| Scaffolds_112 | 664379 | 58.7248469305203 |
| Scaffolds_93 | 659545 | 70.442100250902 |
| Scaffolds_98 | 644476 | 38.2844541409322 |
| Scaffolds_11 | 642189 | 38.6082631480591 |
| Scaffolds_526 | 579093 | 38.3583004515667 |
| Scaffolds_125 | 571100 | 61.7185744776455 |
| Scaffolds_50 | 554491 | 38.8332145642044 |
| Scaffolds_143 | 553638 | 65.0879374796294 |
| Scaffolds_196 | 541223 | 39.4932960800574 |
| Scaffolds_72 | 530323 | 38.3615707349267 |
| Scaffolds_60 | 518110 | 37.6246166396059 |
| Scaffolds_111 | 505903 | 38.9743846708891 |
| Scaffolds_278 | 503841 | 39.8056537453135 |
| Scaffolds_141 | 500204 | 38.2966217918218 |
| Scaffolds_167 | 498008 | 38.5641498271442 |
| Scaffolds_188 | 495100 | 61.0879731866082 |
| Scaffolds_99 | 494368 | 38.4813727382634 |
| Scaffolds_237 | 487549 | 57.8218987234971 |
| Scaffolds_199 | 480716 | 37.6637065702312 |
| Scaffolds_142 | 479561 | 37.719818134983 |
| Scaffolds_104 | 472624 | 60.9409977445168 |
| Scaffolds_84 | 469797 | 62.8040250022461 |
| Scaffolds_45 | 469646 | 39.1730053588212 |
| Scaffolds_7 | 467761 | 37.9036481994164 |
| Scaffolds_827 | 461187 | 47.6058726211049 |
| Scaffolds_781 | 459127 | 38.4839554896336 |
| Scaffolds_71 | 457773 | 37.63530118508 |
| Scaffolds_159 | 441878 | 38.2366253281434 |
| Scaffolds_191 | 441865 | 38.4342410867131 |
| Scaffolds_37 | 439837 | 57.4904250425875 |
| Scaffolds_9 | 438271 | 37.6852503631437 |
| Scaffolds_758 | 432387 | 39.2338843032547 |
| Scaffolds_169 | 429912 | 59.6030559319658 |
| Scaffolds_397 | 423971 | 37.8208046267097 |
| Scaffolds_192 | 421802 | 39.0984740104912 |
| Scaffolds_300 | 419102 | 38.0103745132473 |
| Scaffolds_73 | 417684 | 47.7966036483607 |
| Scaffolds_682 | 403051 | 39.2511654794779 |
| Scaffolds_134 | 396455 | 63.6532736148682 |
| Scaffolds_868 | 392066 | 39.1129135103425 |
| Scaffolds_391 | 380229 | 38.6062682947997 |
| Scaffolds_197 | 375203 | 70.0513347572171 |
| Scaffolds_148 | 370252 | 47.149519785074 |
| Scaffolds_260 | 369669 | 37.6683693366678 |
| Scaffolds_645 | 367440 | 38.2716822155213 |
| Scaffolds_698 | 366300 | 38.6801054867295 |
| Scaffolds_240 | 364443 | 63.7694918602705 |
| Scaffolds_952 | 364443 | 63.7695674244399 |
| Scaffolds_226 | 350636 | 70.8471705775783 |
| Scaffolds_209 | 344756 | 38.4760071609936 |
| Scaffolds_228 | 342101 | 37.9754753657503 |
| Scaffolds_326 | 339259 | 63.3179174735676 |
| Scaffolds_172 | 334139 | 37.9490943658455 |
| Scaffolds_555 | 323962 | 38.0727488239434 |
| Scaffolds_20 | 319779 | 70.1451289572619 |
| Scaffolds_2 | 316177 | 37.9234547518969 |
| Scaffolds_100 | 314385 | 38.4316347485553 |
| Scaffolds_456 | 312949 | 43.5464337867589 |
| Scaffolds_85 | 306253 | 38.2942319309073 |
| Scaffolds_83 | 297684 | 38.7157816761507 |
| Scaffolds_164 | 296665 | 61.9862540406986 |
| Scaffolds_317 | 291447 | 38.0793895329886 |
| Scaffolds_219 | 286800 | 37.3844850555561 |
| Scaffolds_80 | 285074 | 61.3671442001431 |
| Scaffolds_43 | 282425 | 38.0310664348664 |
| Scaffolds_124 | 278038 | 38.8510727206748 |
| Scaffolds_984 | 276524 | 67.228398052986 |
| Scaffolds_31 | 274712 | 38.2937319637783 |
| Scaffolds_225 | 271016 | 63.1747964386272 |
| Scaffolds_234 | 269881 | 63.6366008199776 |
| Scaffolds_66 | 269863 | 37.8760968109641 |
| Scaffolds_170 | 267716 | 39.1257140695343 |
| Scaffolds_208 | 262187 | 38.11948664543 |
| Scaffolds_187 | 260330 | 37.5893090361538 |
| Scaffolds_63 | 259293 | 39.0633145809504 |
| Scaffolds_235 | 249446 | 69.6503394879827 |
| Scaffolds_246 | 247850 | 37.9906320285579 |
| Scaffolds_140 | 246063 | 38.3711620913173 |
| Scaffolds_128 | 243824 | 65.200239190738 |
| Scaffolds_149 | 239730 | 58.2278210466021 |
| Scaffolds_643 | 238000 | 37.4909665383775 |
| Scaffolds_988 | 227786 | 59.5255237325935 |
| Scaffolds_92 | 227001 | 38.1630360564745 |
| Scaffolds_691 | 227001 | 38.1630360564745 |
| Scaffolds_537 | 226299 | 37.6776947516803 |
| Scaffolds_161 | 223170 | 38.8231652639389 |
| Scaffolds_145 | 211448 | 59.9155724628789 |
| Scaffolds_239 | 210643 | 38.3867798395592 |
| Scaffolds_914 | 194764 | 37.9584942084942 |
| Scaffolds_438 | 193899 | 38.2973785591616 |
| Scaffolds_21 | 193151 | 35.9257694093587 |
| Scaffolds_120 | 192530 | 38.2637691396755 |
| Scaffolds_24 | 190405 | 59.7529338440149 |
| Scaffolds_574 | 186373 | 38.4697518443997 |
| Scaffolds_220 | 183373 | 37.6211429664653 |
| Scaffolds_345 | 173949 | 37.8901197449886 |
| Scaffolds_840 | 171420 | 47.1240564221629 |
| Scaffolds_75 | 167010 | 37.9034092950292 |
| Scaffolds_322 | 167010 | 37.9032644360884 |
| Scaffolds_88 | 164384 | 38.1977606356203 |
| Scaffolds_244 | 160833 | 38.492369452321 |
| Scaffolds_656 | 160788 | 37.3992487188417 |
| Scaffolds_633 | 160779 | 37.9404282844366 |
| Scaffolds_1 | 154819 | 38.6187920243378 |
| Scaffolds_105 | 152163 | 39.8446240332432 |
| Scaffolds_981 | 151557 | 68.866250107219 |
| Scaffolds_138 | 148240 | 37.9563898004861 |
| Scaffolds_733 | 146140 | 37.868130063499 |
| Scaffolds_652 | 145142 | 60.3269856142865 |
| Scaffolds_173 | 144802 | 38.5826554148309 |
| Scaffolds_590 | 144802 | 38.5824977210574 |
| Scaffolds_189 | 144240 | 58.8859564604192 |
| Scaffolds_231 | 143470 | 59.4637355897236 |
| Scaffolds_448 | 143470 | 59.4638675142188 |
| Scaffolds_147 | 140779 | 72.1900506792172 |
| Scaffolds_30 | 140549 | 37.7295397465725 |
| Scaffolds_410 | 140549 | 37.7293651414789 |
| Scaffolds_3 | 140191 | 38.1238720906179 |
| Scaffolds_116 | 139317 | 37.5258837178565 |
| Scaffolds_76 | 138639 | 38.9759027054758 |
| Scaffolds_162 | 137344 | 61.1036157582299 |
| Scaffolds_694 | 136172 | 38.8554349582152 |
| Scaffolds_113 | 132439 | 38.7877204868055 |
| Scaffolds_133 | 132104 | 72.2388153551105 |
| Scaffolds_982 | 131873 | 39.5495734597156 |
| Scaffolds_218 | 131451 | 38.5887760538851 |
| Scaffolds_178 | 131254 | 37.9128454579254 |
| Scaffolds_82 | 127224 | 64.5390192629456 |
| Scaffolds_5 | 125264 | 38.0744388826094 |
| Scaffolds_557 | 123366 | 37.9460160492826 |
| Scaffolds_622 | 121638 | 71.8784939164748 |
| Scaffolds_222 | 117403 | 39.1986428346567 |
| Scaffolds_153 | 115109 | 37.7264079643481 |
| Scaffolds_505 | 112811 | 39.0531233102568 |
| Scaffolds_109 | 111827 | 37.4222500714147 |
| Scaffolds_17 | 111538 | 38.1271292809754 |
| Scaffolds_144 | 111156 | 70.760687960688 |
| Scaffolds_118 | 109355 | 36.2820852803204 |
| Scaffolds_56 | 109215 | 61.3512548412793 |
| Scaffolds_132 | 109201 | 47.4797930200937 |
| Scaffolds_212 | 108552 | 38.563992954541 |
| Scaffolds_424 | 108364 | 39.3499806212281 |
| Scaffolds_62 | 107981 | 40.5706851002125 |
| Scaffolds_714 | 107867 | 38.5337770814599 |
| Scaffolds_122 | 107579 | 38.6694821747898 |
| Scaffolds_885 | 105712 | 48.4411063604374 |
| Scaffolds_867 | 105073 | 37.8056092199054 |
| Scaffolds_874 | 104625 | 58.8796390988942 |
| Scaffolds_4 | 104337 | 72.4717424182122 |
| Scaffolds_16 | 104162 | 70.8770578607875 |
| Scaffolds_998 | 103824 | 61.8409646909252 |
| Scaffolds_78 | 102952 | 58.2931041249657 |
| Scaffolds_742 | 101391 | 58.7701563193451 |
| Scaffolds_481 | 99177 | 57.5100069569776 |
| Scaffolds_625 | 98524 | 57.5643471906463 |
| Scaffolds_501 | 97044 | 59.5612480422059 |
| Scaffolds_385 | 96796 | 71.0500216946631 |
| Scaffolds_232 | 96607 | 70.8168164200141 |
| Scaffolds_186 | 96206 | 40.4943539574335 |
| Scaffolds_130 | 91565 | 37.4163181058678 |
| Scaffolds_158 | 89201 | 38.239745300046 |
| Scaffolds_77 | 88722 | 55.9284973626076 |
| Scaffolds_139 | 87221 | 38.4118867729842 |
| Scaffolds_6 | 85106 | 38.5982845729056 |
| Scaffolds_754 | 85006 | 37.8729060794278 |
| Scaffolds_12 | 81392 | 37.1447610834542 |
| Scaffolds_229 | 81042 | 38.7456509321286 |
| Scaffolds_947 | 78867 | 60.6065600791165 |
| Scaffolds_773 | 77169 | 51.6528229516269 |
| Scaffolds_175 | 76459 | 60.6188614100938 |
| Scaffolds_374 | 75402 | 38.2393506975757 |
| Scaffolds_13 | 75098 | 38.8245319698543 |
| Scaffolds_929 | 74792 | 60.4019038960344 |
| Scaffolds_155 | 73930 | 35.7730723132241 |
| Scaffolds_119 | 73423 | 38.2393397524072 |
| Scaffolds_591 | 73357 | 37.9721642879538 |
| Scaffolds_61 | 70971 | 63.9026561804422 |
| Scaffolds_204 | 69890 | 38.8431052737002 |
| Scaffolds_987 | 67846 | 37.5357406042741 |
| Scaffolds_666 | 67645 | 37.8109894008604 |
| Scaffolds_18 | 66613 | 36.8194304757044 |
| Scaffolds_306 | 66591 | 56.971453456069 |
| Scaffolds_32 | 66567 | 38.7862986174233 |
| Scaffolds_25 | 66077 | 38.9037696160772 |
| Scaffolds_948 | 65955 | 60.4226996376427 |
| Scaffolds_248 | 65737 | 63.4774414713002 |
| Scaffolds_110 | 65322 | 39.1508209487962 |
| Scaffolds_152 | 63752 | 38.428961320074 |
| Scaffolds_107 | 63739 | 38.1130256657429 |
| Scaffolds_250 | 62711 | 38.5677947543082 |
| Scaffolds_156 | 62392 | 38.0682198567316 |
| Scaffolds_848 | 61915 | 37.7493176569389 |
| Scaffolds_434 | 61813 | 48.2941033729677 |
| Scaffolds_59 | 61433 | 65.8676932428325 |
| Scaffolds_74 | 60570 | 36.9157507893511 |
| Scaffolds_216 | 60564 | 57.9021233035036 |
| Scaffolds_121 | 60057 | 55.0967897908788 |
| Scaffolds_55 | 59445 | 59.4109038303026 |
| Scaffolds_54 | 59229 | 35.9200959158416 |
| Scaffolds_351 | 58833 | 68.4195319271887 |
| Scaffolds_719 | 58767 | 38.0574792832956 |
| Scaffolds_168 | 57899 | 69.1203329094902 |
| Scaffolds_445 | 57604 | 55.8099569504236 |
| Scaffolds_458 | 57439 | 55.0545758403983 |
| Scaffolds_830 | 57414 | 39.1093075101017 |
| Scaffolds_166 | 57315 | 70.7210926864106 |
| Scaffolds_126 | 57156 | 49.0683229813665 |
| Scaffolds_915 | 56904 | 45.9074295353905 |
| Scaffolds_1017 | 55886 | 48.2500715717149 |
| Scaffolds_51 | 54569 | 39.576347277958 |
| Scaffolds_334 | 54538 | 38.5075174184085 |
| Scaffolds_213 | 54083 | 39.0587953538498 |
| Scaffolds_518 | 54078 | 66.1587959025184 |
| Scaffolds_500 | 53302 | 53.4986680673845 |
| Scaffolds_899 | 53290 | 58.402762140659 |
| Scaffolds_67 | 53212 | 39.2729205211983 |
| Scaffolds_592 | 52955 | 60.5604547085371 |
| Scaffolds_114 | 52897 | 38.7048876056919 |
| Scaffolds_979 | 52406 | 42.6884182407937 |
| Scaffolds_726 | 51999 | 45.4435321039171 |
| Scaffolds_262 | 51745 | 65.9787040794635 |
| Scaffolds_97 | 51654 | 67.739978240857 |
| Scaffolds_876 | 51581 | 55.4601322141015 |
| Scaffolds_227 | 51383 | 39.7821754272945 |
| Scaffolds_38 | 51337 | 38.8323949390936 |
| Scaffolds_47 | 51088 | 38.0157363187975 |
| Scaffolds_353 | 50335 | 32.960527622718 |
| Scaffolds_1024 | 49710 | 47.8396427565676 |
| Scaffolds_684 | 48340 | 39.1713552871091 |
| Scaffolds_459 | 47937 | 38.369594693256 |
| Scaffolds_33 | 47309 | 38.0433113193058 |
| Scaffolds_606 | 47074 | 58.2671311440588 |
| Scaffolds_318 | 46897 | 53.9434127204111 |
| Scaffolds_626 | 46756 | 48.5499807519569 |
| Scaffolds_19 | 46379 | 38.2942026173382 |
| Scaffolds_315 | 46074 | 57.2181952341681 |
| Scaffolds_782 | 45974 | 64.5417373526469 |
| Scaffolds_674 | 45262 | 37.5270622542306 |
| Scaffolds_852 | 44261 | 42.7651643510674 |
| Scaffolds_388 | 44036 | 59.52770208901 |
| Scaffolds_264 | 43945 | 44.465175544381 |
| Scaffolds_973 | 43652 | 68.6205790727506 |
| Scaffolds_261 | 43602 | 38.0309145949913 |
| Scaffolds_223 | 43544 | 71.5751038939201 |
| Scaffolds_478 | 43514 | 45.2180706833954 |
| Scaffolds_1035 | 43280 | 37.6651880602532 |
| Scaffolds_934 | 43225 | 57.7324357461772 |
| Scaffolds_95 | 42883 | 76.321873743466 |
| Scaffolds_48 | 42844 | 38.0904592979836 |
| Scaffolds_854 | 42364 | 67.1599867818534 |
| Scaffolds_27 | 41945 | 36.5348399246704 |
| Scaffolds_301 | 41490 | 48.4335084590543 |
| Scaffolds_900 | 41337 | 38.4654459253526 |
| Scaffolds_328 | 40729 | 47.5560356467729 |
| Scaffolds_717 | 40682 | 35.5643496214728 |
| Scaffolds_179 | 40215 | 51.8402376629807 |
| Scaffolds_426 | 40071 | 48.613849033063 |
| Scaffolds_895 | 39947 | 40.6948511927111 |
| Scaffolds_290 | 39613 | 55.0271361857882 |
| Scaffolds_256 | 39310 | 38.9041514041514 |
| Scaffolds_154 | 39251 | 36.0038023817697 |
| Scaffolds_26 | 38818 | 38.4575756014631 |
| Scaffolds_414 | 38645 | 48.7929830008538 |
| Scaffolds_406 | 38417 | 41.6981181186392 |
| Scaffolds_34 | 38181 | 72.4865784994108 |
| Scaffolds_917 | 37883 | 66.2654736453137 |
| Scaffolds_470 | 37583 | 66.2312092590129 |
| Scaffolds_215 | 37304 | 37.669065748757 |
| Scaffolds_79 | 36965 | 64.1862100678947 |
| Scaffolds_761 | 36796 | 48.0190217391304 |
| Scaffolds_785 | 36690 | 35.5881391038919 |
| Scaffolds_198 | 36572 | 37.3712431212078 |
| Scaffolds_108 | 36458 | 47.2940211107381 |
| Scaffolds_90 | 36408 | 38.0225212853612 |
| Scaffolds_243 | 36383 | 36.7988567345481 |
| Scaffolds_324 | 35664 | 37.9717377894914 |
| Scaffolds_245 | 35486 | 44.6790006325111 |
| Scaffolds_702 | 35142 | 43.310380150239 |
| Scaffolds_151 | 35084 | 39.4288320127686 |
| Scaffolds_194 | 35062 | 37.910447761194 |
| Scaffolds_467 | 34799 | 35.7670315777376 |
| Scaffolds_283 | 34797 | 58.3833797879371 |
| Scaffolds_844 | 34432 | 60.3902892321989 |
| Scaffolds_760 | 34065 | 69.3504358801256 |
| Scaffolds_330 | 33737 | 39.0089208974245 |
| Scaffolds_428 | 33378 | 37.2266490923252 |
| Scaffolds_490 | 33305 | 42.637125101324 |
| Scaffolds_460 | 33105 | 37.6536684084937 |
| Scaffolds_180 | 32862 | 37.8932635550417 |
| Scaffolds_880 | 32862 | 37.8925267770204 |
| Scaffolds_182 | 32859 | 43.4684982168802 |
| Scaffolds_49 | 32796 | 56.681941457028 |
| Scaffolds_683 | 32764 | 56.9767441860465 |
| Scaffolds_65 | 32632 | 57.1982660278348 |
| Scaffolds_597 | 32313 | 56.9947705541975 |
| Scaffolds_743 | 32275 | 67.8273800303603 |
| Scaffolds_1013 | 32210 | 58.4347448155967 |
| Scaffolds_64 | 32001 | 36.1877793319914 |
| Scaffolds_255 | 31225 | 37.9006692497358 |
| Scaffolds_529 | 31225 | 37.9006692497358 |
| Scaffolds_832 | 31210 | 34.6254885628244 |
| Scaffolds_941 | 30778 | 37.7745289148798 |
| Scaffolds_131 | 30539 | 36.2538061094195 |
| Scaffolds_863 | 30375 | 55.2457451361227 |
| Scaffolds_247 | 29941 | 22.7005737955387 |
| Scaffolds_342 | 29744 | 36.1805835686433 |
| Scaffolds_313 | 29438 | 58.4233408056518 |
| Scaffolds_909 | 29325 | 55.7571004807528 |
| Scaffolds_279 | 29159 | 57.7306861434009 |
| Scaffolds_544 | 29049 | 47.9606236877431 |
| Scaffolds_267 | 29042 | 35.5505060937823 |
| Scaffolds_346 | 28888 | 36.2166839736933 |
| Scaffolds_89 | 28527 | 38.245294261979 |
| Scaffolds_86 | 28469 | 36.7857268289256 |
| Scaffolds_302 | 28400 | 54.8514293761442 |
| Scaffolds_361 | 28185 | 59.2358721487105 |
| Scaffolds_1016 | 27895 | 47.3923796551848 |
| Scaffolds_1042 | 27880 | 36.2860421747239 |
| Scaffolds_685 | 27793 | 57.1428571428571 |
| Scaffolds_786 | 27621 | 35.3122171945701 |
| Scaffolds_339 | 27549 | 42.1400312148379 |
| Scaffolds_15 | 27239 | 36.6477272727273 |
| Scaffolds_276 | 27174 | 39.4392110685899 |
| Scaffolds_207 | 27040 | 35.5311504858333 |
| Scaffolds_545 | 26782 | 36.1820489844683 |
| Scaffolds_41 | 26730 | 42.8861887011667 |
| Scaffolds_962 | 26670 | 53.5727674889405 |
| Scaffolds_986 | 26647 | 54.6095831300889 |
| Scaffolds_549 | 26531 | 36.6007160354249 |
| Scaffolds_706 | 26515 | 40.3899091217618 |
| Scaffolds_440 | 26510 | 36.3619220034699 |
| Scaffolds_745 | 26481 | 42.7373985274684 |
| Scaffolds_320 | 26362 | 58.5033755594326 |
| Scaffolds_834 | 26287 | 35.0728386139744 |
| Scaffolds_297 | 26271 | 41.1379638439581 |
| Scaffolds_493 | 26186 | 63.4718191538109 |
| Scaffolds_68 | 25947 | 37.1978953000981 |
| Scaffolds_514 | 25927 | 49.4388955304462 |
| Scaffolds_847 | 25914 | 61.6723259762309 |
| Scaffolds_190 | 25893 | 69.4698099142751 |
| Scaffolds_185 | 25531 | 37.3577786948063 |
| Scaffolds_233 | 25436 | 59.2029673858959 |
| Scaffolds_238 | 25380 | 37.7481878348566 |
| Scaffolds_303 | 25380 | 37.7472224411 |
| Scaffolds_39 | 25327 | 56.913655228262 |
| Scaffolds_338 | 25252 | 38.3314855875831 |
| Scaffolds_419 | 25238 | 56.3386419459631 |
| Scaffolds_277 | 25102 | 58.3963992671075 |
| Scaffolds_820 | 25023 | 58.688616294402 |
| Scaffolds_506 | 24907 | 67.9152113693846 |
| Scaffolds_205 | 24886 | 37.6657292085175 |
| Scaffolds_280 | 24758 | 57.2893950407883 |
| Scaffolds_367 | 24734 | 37.5454765947126 |
| Scaffolds_513 | 24573 | 70.858933148879 |
| Scaffolds_766 | 24543 | 37.6451415766959 |
| Scaffolds_201 | 24441 | 38.3211862276954 |
| Scaffolds_510 | 24408 | 48.7342290676716 |
| Scaffolds_344 | 24405 | 53.8653775246835 |
| Scaffolds_422 | 24167 | 38.6983325747859 |
| Scaffolds_91 | 24122 | 37.5652822680925 |
| Scaffolds_340 | 24101 | 57.9838207840697 |
| Scaffolds_556 | 23993 | 42.6963950823088 |
| Scaffolds_29 | 23980 | 36.899051868377 |
| Scaffolds_752 | 23705 | 42.9288455860644 |
| Scaffolds_855 | 23335 | 59.0899353014268 |
| Scaffolds_916 | 23304 | 38.1757336536811 |
| Scaffolds_437 | 23164 | 56.9147099447514 |
| Scaffolds_935 | 23089 | 48.7247217771619 |
| Scaffolds_650 | 23083 | 52.7312107429067 |
| Scaffolds_171 | 23012 | 38.3255126868266 |
| Scaffolds_357 | 23012 | 38.3255126868266 |
| Scaffolds_568 | 22976 | 62.2541340295909 |
| Scaffolds_668 | 22804 | 57.1159242371098 |
| Scaffolds_731 | 22677 | 38.1057365845055 |
| Scaffolds_252 | 22306 | 38.2294935006723 |
| Scaffolds_776 | 22306 | 38.2294935006723 |
| Scaffolds_288 | 22292 | 48.0400071761751 |
| Scaffolds_771 | 22208 | 47.3257698541329 |
| Scaffolds_417 | 22129 | 58.0129218813536 |
| Scaffolds_747 | 22085 | 47.9061976549414 |
| Scaffolds_362 | 22077 | 45.3059191159821 |
| Scaffolds_408 | 22059 | 43.4890993971808 |
| Scaffolds_291 | 21968 | 42.8818496267977 |
| Scaffolds_823 | 21959 | 63.7116969448618 |
| Scaffolds_808 | 21934 | 36.5165466314158 |
| Scaffolds_757 | 21815 | 44.8416517713919 |
| Scaffolds_787 | 21765 | 58.5473423071622 |
| Scaffolds_35 | 21636 | 70.3234750462107 |
| Scaffolds_347 | 21628 | 37.0192307692308 |
| Scaffolds_485 | 21606 | 36.682091624248 |
| Scaffolds_462 | 21400 | 58.3769388899271 |
| Scaffolds_965 | 21368 | 53.0787946846341 |
| Scaffolds_980 | 21367 | 26.7090917598615 |
| Scaffolds_271 | 21329 | 48.0335630244223 |
| Scaffolds_583 | 21175 | 59.3040275744842 |
| Scaffolds_901 | 21093 | 41.8590320898706 |
| Scaffolds_377 | 21089 | 33.3949651543166 |
| Scaffolds_696 | 21028 | 58.6772537086345 |
| Scaffolds_535 | 21026 | 40.190204469805 |
| Scaffolds_532 | 20981 | 48.9063616869192 |
| Scaffolds_693 | 20962 | 69.2454450062005 |
| Scaffolds_394 | 20911 | 42.648816638776 |
| Scaffolds_833 | 20848 | 48.9545367350854 |
| Scaffolds_519 | 20795 | 38.0258666282033 |
| Scaffolds_887 | 20713 | 44.2824733310808 |
| Scaffolds_921 | 20530 | 51.694750170449 |
| Scaffolds_669 | 20521 | 58.9866017052375 |
| Scaffolds_858 | 20521 | 38.2898903775883 |
| Scaffolds_588 | 20469 | 31.5146778684121 |
| Scaffolds_646 | 20468 | 71.73700664322 |
| Scaffolds_603 | 20262 | 40.4766604164611 |
| Scaffolds_775 | 20227 | 57.1548613513914 |
| Scaffolds_558 | 19941 | 23.9107545750815 |
| Scaffolds_727 | 19887 | 58.2172841988839 |
| Scaffolds_1037 | 19849 | 37.4855185614265 |
| Scaffolds_632 | 19767 | 37.9040008092661 |
| Scaffolds_736 | 19610 | 58.8559192413582 |
| Scaffolds_751 | 19523 | 38.72074563425 |
| Scaffolds_580 | 19393 | 37.7584162499356 |
| Scaffolds_503 | 19381 | 55.2953314418365 |
| Scaffolds_849 | 19320 | 39.2517077209687 |
| Scaffolds_411 | 19205 | 38.7786974855536 |
| Scaffolds_364 | 19132 | 56.5844481605351 |
| Scaffolds_755 | 19037 | 38.280552491991 |
| Scaffolds_601 | 18964 | 47.6328553353016 |
| Scaffolds_695 | 18913 | 46.926045356029 |
| Scaffolds_951 | 18860 | 35.437871077184 |
| Scaffolds_681 | 18840 | 59.8864360008491 |
| Scaffolds_383 | 18771 | 62.982689747004 |
| Scaffolds_671 | 18721 | 43.5567423230975 |
| Scaffolds_657 | 18584 | 69.1629007962126 |
| Scaffolds_354 | 18516 | 40.1825251107031 |
| Scaffolds_36 | 18480 | 73.2958234148453 |
| Scaffolds_443 | 18472 | 46.6334704481489 |
| Scaffolds_843 | 18372 | 42.1745755333043 |
| Scaffolds_919 | 18351 | 35.9738490874421 |
| Scaffolds_711 | 18254 | 35.9568408368934 |
| Scaffolds_531 | 18252 | 48.1102103418054 |
| Scaffolds_658 | 18196 | 42.032967032967 |
| Scaffolds_879 | 18191 | 52.3110744710085 |
| Scaffolds_534 | 18076 | 43.3185840707965 |
| Scaffolds_664 | 18032 | 41.9827012641384 |
| Scaffolds_349 | 17916 | 47.7006362317223 |
| Scaffolds_117 | 17872 | 36.9792433387311 |
| Scaffolds_157 | 17870 | 39.1294617880721 |
| Scaffolds_495 | 17650 | 37.9687322986292 |
| Scaffolds_627 | 17633 | 62.0513692804899 |
| Scaffolds_564 | 17627 | 48.306959332993 |
| Scaffolds_287 | 17442 | 36.5642554167144 |
| Scaffolds_477 | 17426 | 52.6965002868617 |
| Scaffolds_974 | 17380 | 58.8989875747814 |
| Scaffolds_911 | 17379 | 66.7606281999655 |
| Scaffolds_819 | 17309 | 49.298215214001 |
| Scaffolds_450 | 17245 | 46.5824105745261 |
| Scaffolds_686 | 17220 | 71.8590339061774 |
| Scaffolds_380 | 17123 | 40.0362001518071 |
| Scaffolds_373 | 17034 | 39.0949641976758 |
| Scaffolds_737 | 17028 | 43.8637698179683 |
| Scaffolds_176 | 16961 | 66.5546713822576 |
| Scaffolds_689 | 16941 | 56.4060194747713 |
| Scaffolds_577 | 16911 | 47.0469997044044 |
| Scaffolds_273 | 16900 | 53.2122574538571 |
| Scaffolds_792 | 16881 | 59.3081798258603 |
| Scaffolds_94 | 16839 | 55.6314195808348 |
| Scaffolds_593 | 16750 | 56.4641279694401 |
| Scaffolds_309 | 16705 | 57.7293674067868 |
| Scaffolds_992 | 16653 | 36.0449060455064 |
| Scaffolds_497 | 16523 | 61.4993646759847 |
| Scaffolds_546 | 16522 | 68.1493585088356 |
| Scaffolds_115 | 16482 | 21.8376125823363 |
| Scaffolds_217 | 16403 | 36.9893428063943 |
| Scaffolds_366 | 16383 | 43.9982913284921 |
| Scaffolds_605 | 16375 | 58.25141950058 |
| Scaffolds_985 | 16279 | 34.1337591352945 |
| Scaffolds_619 | 16198 | 58.3641975308642 |
| Scaffolds_230 | 16115 | 54.2403781228899 |
| Scaffolds_898 | 16099 | 56.7347699186487 |
| Scaffolds_765 | 16052 | 38.6646736422521 |
| Scaffolds_371 | 16051 | 70.9000311429461 |
| Scaffolds_660 | 16019 | 46.9075703675966 |
| Scaffolds_489 | 16005 | 58.2234992816541 |
| Scaffolds_471 | 15991 | 58.6245701781807 |
| Scaffolds_824 | 15944 | 59.8883872585904 |
| Scaffolds_457 | 15846 | 58.2397476340694 |
| Scaffolds_553 | 15826 | 60.7454200884397 |
| Scaffolds_304 | 15817 | 37.7283357562733 |
| Scaffolds_896 | 15726 | 48.5187539732994 |
| Scaffolds_614 | 15601 | 33.2906119833387 |
| Scaffolds_735 | 15583 | 34.9714505677808 |
| Scaffolds_648 | 15412 | 45.9782044628957 |
| Scaffolds_883 | 15387 | 37.4569553635241 |
| Scaffolds_990 | 15371 | 63.0178861788618 |
| Scaffolds_953 | 15327 | 34.263909725393 |
| Scaffolds_305 | 15311 | 57.7407770159974 |
| Scaffolds_400 | 15250 | 59.5253703946506 |
| Scaffolds_451 | 15249 | 36.1174850849013 |
| Scaffolds_651 | 15206 | 49.2766964755392 |
| Scaffolds_242 | 15137 | 76.5008916187834 |
| Scaffolds_407 | 15137 | 76.5008916187834 |
| Scaffolds_1004 | 14993 | 39.8413015936521 |
| Scaffolds_926 | 14987 | 50.8438396371156 |
| Scaffolds_778 | 14977 | 59.7623656631734 |
| Scaffolds_359 | 14976 | 57.3898531375167 |
| Scaffolds_418 | 14954 | 43.1073672950929 |
| Scaffolds_728 | 14869 | 63.6657029516574 |
| Scaffolds_913 | 14852 | 43.3360258481422 |
| Scaffolds_956 | 14830 | 47.0334412081985 |
| Scaffolds_127 | 14829 | 35.4729170117608 |
| Scaffolds_1009 | 14825 | 21.4242362937487 |
| Scaffolds_732 | 14783 | 42.4494488401975 |
| Scaffolds_594 | 14781 | 48.0825160635779 |
| Scaffolds_969 | 14665 | 35.9124684709251 |
| Scaffolds_963 | 14657 | 36.127975987448 |
| Scaffolds_214 | 14596 | 35.6763817454486 |
| Scaffolds_559 | 14552 | 42.9032701291564 |
| Scaffolds_402 | 14527 | 38.2974330741174 |
| Scaffolds_386 | 14517 | 55.3336547069761 |
| Scaffolds_875 | 14476 | 41.2016574585635 |
| Scaffolds_274 | 14464 | 37.0472767486868 |
| Scaffolds_763 | 14419 | 34.8609859252583 |
| Scaffolds_807 | 14331 | 39.2396232996163 |
| Scaffolds_515 | 14254 | 57.7991303128068 |
| Scaffolds_565 | 14227 | 36.7507553931558 |
| Scaffolds_630 | 14210 | 59.4638333802421 |
| Scaffolds_241 | 14208 | 56.1237928007024 |
| Scaffolds_616 | 14208 | 69.4862772695285 |
| Scaffolds_46 | 14207 | 41.0738160579833 |
| Scaffolds_925 | 14182 | 35.1825743690963 |
| Scaffolds_1029 | 14149 | 50.7524906380273 |
| Scaffolds_908 | 14125 | 69.2688796093142 |
| Scaffolds_634 | 14083 | 59.9062965855044 |
| Scaffolds_541 | 14040 | 38.5075477072059 |
| Scaffolds_803 | 14029 | 59.3885840518777 |
| Scaffolds_1007 | 13975 | 47.4851541818702 |
| Scaffolds_203 | 13840 | 35.2354810748339 |
| Scaffolds_379 | 13830 | 55.3129969640017 |
| Scaffolds_352 | 13828 | 37.2831116252169 |
| Scaffolds_413 | 13770 | 35.1168868883404 |
| Scaffolds_759 | 13756 | 37.4273255813954 |
| Scaffolds_809 | 13756 | 57.703488372093 |
| Scaffolds_665 | 13749 | 46.6443685014179 |
| Scaffolds_523 | 13736 | 70.8588064046579 |
| Scaffolds_586 | 13729 | 54.0231559018423 |
| Scaffolds_642 | 13724 | 56.2791375291375 |
| Scaffolds_599 | 13710 | 61.3679451655243 |
| Scaffolds_994 | 13653 | 37.0139855019404 |
| Scaffolds_931 | 13604 | 60.7289829512052 |
| Scaffolds_395 | 13601 | 39.1604793060354 |
| Scaffolds_163 | 13592 | 65.0338334804354 |
| Scaffolds_551 | 13576 | 45.4565537555228 |
| Scaffolds_511 | 13564 | 41.6199882075472 |
| Scaffolds_569 | 13548 | 47.6387249114522 |
| Scaffolds_722 | 13361 | 35.9221848110737 |
| Scaffolds_562 | 13346 | 61.8801498127341 |
| Scaffolds_598 | 13338 | 42.2200569629741 |
| Scaffolds_1025 | 13292 | 67.2157039711191 |
| Scaffolds_610 | 13197 | 26.3086129838649 |
| Scaffolds_508 | 13129 | 63.6640523871164 |
| Scaffolds_533 | 12939 | 37.8505756007108 |
| Scaffolds_730 | 12936 | 38.0602782071097 |
| Scaffolds_850 | 12907 | 36.263651150182 |
| Scaffolds_764 | 12867 | 48.6597777950431 |
| Scaffolds_368 | 12748 | 38.5664993726474 |
| Scaffolds_447 | 12717 | 47.8185677226633 |
| Scaffolds_949 | 12696 | 55.3385826771653 |
| Scaffolds_292 | 12692 | 24.9448645242596 |
| Scaffolds_810 | 12631 | 58.2904629996043 |
| Scaffolds_265 | 12575 | 39.0889577867875 |
| Scaffolds_888 | 12563 | 47.1552478714092 |
| Scaffolds_976 | 12528 | 35.8123204596234 |
| Scaffolds_955 | 12521 | 44.3193612774451 |
| Scaffolds_741 | 12481 | 21.057268722467 |
| Scaffolds_636 | 12453 | 59.918118327045 |
| Scaffolds_480 | 12436 | 67.3580961569384 |
| Scaffolds_183 | 12330 | 58.454647256439 |
| Scaffolds_910 | 12181 | 57.3574066475174 |
| Scaffolds_748 | 12062 | 62.5476545665506 |
| Scaffolds_491 | 12051 | 37.2604330871982 |
| Scaffolds_57 | 11998 | 33.2361273121146 |
| Scaffolds_1020 | 11958 | 59.7475338572145 |
| Scaffolds_624 | 11857 | 47.1714020740241 |
| Scaffolds_463 | 11828 | 67.2836375929682 |
| Scaffolds_699 | 11805 | 58.5570327716149 |
| Scaffolds_928 | 11751 | 41.8885580603998 |
| Scaffolds_663 | 11730 | 54.6275779785239 |
| Scaffolds_521 | 11686 | 39.5295124037639 |
| Scaffolds_677 | 11663 | 41.7759492585926 |
| Scaffolds_482 | 11518 | 37.6931088352717 |
| Scaffolds_433 | 11480 | 56.1389759665622 |
| Scaffolds_893 | 11370 | 39.0715667311412 |
| Scaffolds_360 | 11331 | 62.4614027348919 |
| Scaffolds_517 | 11313 | 58.4165414862596 |
| Scaffolds_412 | 11283 | 36.3249756356871 |
| Scaffolds_975 | 11238 | 33.8551859099804 |
| Scaffolds_316 | 11185 | 33.5418714809188 |
| Scaffolds_611 | 11136 | 36.5798922800718 |
| Scaffolds_509 | 11131 | 46.3044454422991 |
| Scaffolds_582 | 11122 | 60.7855473665288 |
| Scaffolds_886 | 11095 | 41.1568609784665 |
| Scaffolds_708 | 11057 | 59.7775969623 |
| Scaffolds_390 | 11037 | 46.6443256951363 |
| Scaffolds_842 | 11017 | 36.8569095363397 |
| Scaffolds_319 | 11013 | 36.0533720613597 |
| Scaffolds_825 | 11013 | 56.0769719524371 |
| Scaffolds_525 | 10985 | 57.0674433421316 |
| Scaffolds_87 | 10970 | 35.9486057955167 |
| Scaffolds_449 | 10966 | 26.1166818596171 |
| Scaffolds_839 | 10961 | 36.5709074327405 |
| Scaffolds_930 | 10938 | 38.2745384755986 |
| Scaffolds_268 | 10929 | 58.9499679868289 |
| Scaffolds_784 | 10899 | 36.8705860772265 |
| Scaffolds_299 | 10864 | 34.1645196908355 |
| Scaffolds_479 | 10828 | 59.3611521418021 |
| Scaffolds_1005 | 10773 | 35.6778324208964 |
| Scaffolds_972 | 10761 | 38.6065954482118 |
| Scaffolds_966 | 10744 | 34.9925567547451 |
| Scaffolds_971 | 10742 | 37.8652521868602 |
| Scaffolds_615 | 10685 | 49.6304612218168 |
| Scaffolds_692 | 10615 | 61.2204539033807 |
| Scaffolds_355 | 10556 | 41.0890151515151 |
| Scaffolds_799 | 10554 | 42.0344762265581 |
| Scaffolds_476 | 10513 | 55.3960254825521 |
| Scaffolds_811 | 10486 | 58.3222116301239 |
| Scaffolds_851 | 10477 | 37.2578952390039 |
| Scaffolds_673 | 10474 | 56.4611567092957 |
| Scaffolds_983 | 10416 | 34.7120921305182 |
| Scaffolds_293 | 10397 | 42.5343716950293 |
| Scaffolds_573 | 10335 | 68.7687397233775 |
| Scaffolds_640 | 10133 | 59.1200552431686 |
| Scaffolds_581 | 10113 | 60.5614312543244 |
| Scaffolds_1031 | 10068 | 55.8876092136616 |
| Scaffolds_454 | 10043 | 54.1156564148502 |
| Scaffolds_738 | 10014 | 35.436214813336 |
| Scaffolds_923 | 10006 | 69.3606393606394 |
| Scaffolds_836 | 9925 | 61.4361969986907 |
| Scaffolds_499 | 9874 | 67.1998380238915 |
| Scaffolds_403 | 9836 | 54.4105691056911 |
| Scaffolds_1000 | 9834 | 60.9981703598292 |
| Scaffolds_922 | 9830 | 51.6168395363026 |
| Scaffolds_637 | 9820 | 56.3517915309446 |
| Scaffolds_356 | 9802 | 37.4770548643688 |
| Scaffolds_884 | 9783 | 45.3663022376622 |
| Scaffolds_817 | 9762 | 61.1714110178169 |
| Scaffolds_589 | 9727 | 39.029904429144 |
| Scaffolds_802 | 9636 | 56.0373443983403 |
| Scaffolds_1012 | 9635 | 58.1180620396307 |
| Scaffolds_670 | 9626 | 32.3883696780893 |
| Scaffolds_939 | 9575 | 57.6469360058461 |
| Scaffolds_889 | 9560 | 65.7360936846508 |
| Scaffolds_566 | 9558 | 65.4465592972182 |
| Scaffolds_749 | 9550 | 59.4096713418463 |
| Scaffolds_275 | 9507 | 70.6234885921564 |
| Scaffolds_873 | 9462 | 32.2945277836467 |
| Scaffolds_504 | 9447 | 58.0256057560047 |
| Scaffolds_734 | 9368 | 34.816474605207 |
| Scaffolds_560 | 9362 | 55.2743967542174 |
| Scaffolds_712 | 9339 | 62.9776303114631 |
| Scaffolds_442 | 9333 | 56.9561957802292 |
| Scaffolds_392 | 9265 | 40.6624231308663 |
| Scaffolds_762 | 9230 | 36.6904916612519 |
| Scaffolds_991 | 9209 | 36.144578313253 |
| Scaffolds_365 | 9146 | 44.8196721311475 |
| Scaffolds_933 | 9125 | 36.2909409573885 |
| Scaffolds_425 | 9084 | 49.449823943662 |
| Scaffolds_723 | 9084 | 59.3419894366197 |
| Scaffolds_286 | 9074 | 60.2775941837409 |
| Scaffolds_554 | 9055 | 53.1957169665526 |
| Scaffolds_1011 | 9039 | 59.5267057392458 |
| Scaffolds_335 | 9029 | 41.4369533931141 |
| Scaffolds_769 | 8968 | 35.8448506464556 |
| Scaffolds_404 | 8946 | 40.4134078212291 |
| Scaffolds_890 | 8944 | 71.4573088958426 |
| Scaffolds_435 | 8923 | 61.0955528172958 |
| Scaffolds_575 | 8923 | 53.8702811694858 |
| Scaffolds_687 | 8902 | 66.7752077251291 |
| Scaffolds_436 | 8872 | 32.6611086074808 |
| Scaffolds_608 | 8848 | 58.2354270221419 |
| Scaffolds_835 | 8846 | 46.9152542372881 |
| Scaffolds_845 | 8813 | 46.8753544289441 |
| Scaffolds_272 | 8808 | 36.9723104857013 |
| Scaffolds_538 | 8781 | 43.8019351166762 |
| Scaffolds_801 | 8680 | 60.5711653615845 |
| Scaffolds_902 | 8678 | 38.2745911080396 |
| Scaffolds_576 | 8654 | 49.0644490644491 |
| Scaffolds_860 | 8654 | 59.1245091245091 |
| Scaffolds_446 | 8587 | 38.4588522872774 |
| Scaffolds_563 | 8576 | 63.3682983682984 |
| Scaffolds_797 | 8498 | 37.3206304398965 |
| Scaffolds_269 | 8497 | 69.5212327961416 |
| Scaffolds_718 | 8466 | 35.6434474616293 |
| Scaffolds_829 | 8349 | 37.7469172752305 |
| Scaffolds_536 | 8294 | 44.4805977343938 |
| Scaffolds_789 | 8280 | 73.3099951714148 |
| Scaffolds_332 | 8246 | 37.2 |
| Scaffolds_653 | 8230 | 46.7087685207675 |
| Scaffolds_427 | 8220 | 59.7641050583658 |
| Scaffolds_912 | 8220 | 36.8920233463035 |
| Scaffolds_415 | 8190 | 67.280937271174 |
| Scaffolds_672 | 8171 | 34.5932721712538 |
| Scaffolds_486 | 8168 | 49.0210474791973 |
| Scaffolds_767 | 8148 | 55.5569185475957 |
| Scaffolds_993 | 8134 | 57.5694273777341 |
| Scaffolds_870 | 8087 | 38.8579903596589 |
| Scaffolds_542 | 8078 | 56.1494679534769 |
| Scaffolds_298 | 8070 | 60.6267029972752 |
| Scaffolds_572 | 8037 | 46.1634125108817 |
| Scaffolds_826 | 8024 | 54.5092177379173 |
| Scaffolds_805 | 8012 | 61.6641716566866 |
| Scaffolds_469 | 7999 | 42.3466200174934 |
| Scaffolds_788 | 7998 | 59.9850037490627 |
| Scaffolds_1008 | 7990 | 34.4008006004503 |
| Scaffolds_724 | 7977 | 67.1720335797519 |
| Scaffolds_370 | 7976 | 37.3934837092732 |
| Scaffolds_1036 | 7965 | 35.5753544986824 |
| Scaffolds_794 | 7964 | 33.2078313253012 |
| Scaffolds_996 | 7962 | 36.5302535777052 |
| Scaffolds_602 | 7955 | 56.9795200402061 |
| Scaffolds_944 | 7955 | 60.6985802236462 |
| Scaffolds_444 | 7954 | 62.8298567479266 |
| Scaffolds_957 | 7945 | 53.5790665492515 |
| Scaffolds_484 | 7935 | 37.3346769114498 |
| Scaffolds_416 | 7931 | 51.8966603654694 |
| Scaffolds_796 | 7898 | 64.7684130599848 |
| Scaffolds_864 | 7896 | 70.4050632911392 |
| Scaffolds_628 | 7849 | 59.6205271870623 |
| Scaffolds_977 | 7787 | 38.9680400462072 |
| Scaffolds_550 | 7786 | 39.9229781771502 |
| Scaffolds_713 | 7756 | 48.7628865979381 |
| Scaffolds_331 | 7754 | 46.2490332559938 |
| Scaffolds_853 | 7738 | 42.1338155515371 |
| Scaffolds_488 | 7710 | 66.0487425460202 |
| Scaffolds_768 | 7709 | 35.4596136393103 |
| Scaffolds_570 | 7672 | 34.4580510682647 |
| Scaffolds_578 | 7662 | 58.1659274719541 |
| Scaffolds_701 | 7656 | 38.8381201044386 |
| Scaffolds_547 | 7630 | 37.5818705789887 |
| Scaffolds_861 | 7625 | 67.4793550924105 |
| Scaffolds_750 | 7614 | 35.416119716461 |
| Scaffolds_790 | 7554 | 61.0082032283673 |
| Scaffolds_1014 | 7549 | 69.0454124189064 |
| Scaffolds_487 | 7488 | 59.1430859583556 |
| Scaffolds_905 | 7474 | 49.3044408774746 |
| Scaffolds_70 | 7470 | 38.4131656408884 |
| Scaffolds_655 | 7460 | 65.058949624866 |
| Scaffolds_333 | 7456 | 47.8954423592493 |
| Scaffolds_680 | 7373 | 24.1426053951471 |
| Scaffolds_310 | 7366 | 58.5074626865672 |
| Scaffolds_964 | 7326 | 37.7456331877729 |
| Scaffolds_257 | 7322 | 34.015834015834 |
| Scaffolds_372 | 7322 | 34.0114691425451 |
| Scaffolds_102 | 7254 | 62.2794928335171 |
| Scaffolds_206 | 7231 | 38.2031789910159 |
| Scaffolds_539 | 7218 | 57.7125450013847 |
| Scaffolds_989 | 7202 | 36.8581737441021 |
| Scaffolds_329 | 7193 | 53.6612477421148 |
| Scaffolds_1033 | 7141 | 64.0587823652904 |
| Scaffolds_1030 | 7122 | 59.865282065675 |
| Scaffolds_920 | 7116 | 34.5084269662921 |
| Scaffolds_327 | 7074 | 41.8762362249223 |
| Scaffolds_705 | 7061 | 36.0509554140127 |
| Scaffolds_903 | 7053 | 44.8349156865524 |
| Scaffolds_584 | 7027 | 44.1900156450007 |
| Scaffolds_882 | 7027 | 38.515147205234 |
| Scaffolds_716 | 7016 | 72.7777777777778 |
| Scaffolds_421 | 7015 | 65.2657073657216 |
| Scaffolds_567 | 7001 | 37.4732334047109 |
| Scaffolds_924 | 6965 | 33.5342229875161 |
| Scaffolds_548 | 6912 | 60.1792943898207 |
| Scaffolds_1006 | 6908 | 75.853587962963 |
| Scaffolds_507 | 6896 | 48.9130434782609 |
| Scaffolds_881 | 6887 | 60.1073864460891 |
| Scaffolds_1023 | 6838 | 33.8351359251681 |
| Scaffolds_323 | 6817 | 35.2294384987539 |
| Scaffolds_697 | 6817 | 38.8945902360358 |
| Scaffolds_314 | 6792 | 63.3460859329017 |
| Scaffolds_378 | 6790 | 48.0718280836032 |
| Scaffolds_623 | 6790 | 38.519281719164 |
| Scaffolds_635 | 6749 | 61.8836072856508 |
| Scaffolds_950 | 6721 | 36.6988847583643 |
| Scaffolds_1038 | 6720 | 37.5371802498513 |
| Scaffolds_800 | 6706 | 76.1102831594635 |
| Scaffolds_772 | 6699 | 62.4645681038341 |
| Scaffolds_779 | 6676 | 39.940119760479 |
| Scaffolds_936 | 6639 | 54.930001505344 |
| Scaffolds_904 | 6633 | 41.1631761337954 |
| Scaffolds_384 | 6595 | 56.8267919381724 |
| Scaffolds_527 | 6558 | 36.6351722035965 |
| Scaffolds_918 | 6557 | 31.9615912208505 |
| Scaffolds_452 | 6543 | 57.5072552314037 |
| Scaffolds_667 | 6542 | 44.8212648945921 |
| Scaffolds_321 | 6540 | 40.7243276283619 |
| Scaffolds_846 | 6518 | 36.5685372585097 |
| Scaffolds_894 | 6496 | 63.3076923076923 |
| Scaffolds_455 | 6455 | 66.8989007586314 |
| Scaffolds_621 | 6410 | 32.8344246959775 |
| Scaffolds_700 | 6365 | 50.196263149631 |
| Scaffolds_494 | 6356 | 40.3616352201258 |
| Scaffolds_739 | 6353 | 62.5924178071417 |
| Scaffolds_282 | 6352 | 66.9918187539333 |
| Scaffolds_812 | 6343 | 37.8919174413109 |
| Scaffolds_381 | 6315 | 34.8947618294034 |
| Scaffolds_1032 | 6315 | 35.7018515587909 |
| Scaffolds_296 | 6296 | 47.4444444444444 |
| Scaffolds_1010 | 6278 | 37.2652021649156 |
| Scaffolds_461 | 6246 | 32.496 |
| Scaffolds_496 | 6203 | 37.8765909457065 |
| Scaffolds_828 | 6186 | 33.3279483037157 |
| Scaffolds_307 | 6180 | 33.1824062095731 |
| Scaffolds_869 | 6168 | 59.5430978613091 |
| Scaffolds_595 | 6162 | 57.5737917612715 |
| Scaffolds_720 | 6159 | 36.6217751095246 |
| Scaffolds_524 | 6151 | 38.4402924451665 |
| Scaffolds_399 | 6146 | 37.8048780487805 |
| Scaffolds_358 | 6143 | 38.4577842850171 |
| Scaffolds_1019 | 6143 | 34.3419554254108 |
| Scaffolds_475 | 6141 | 64.8543057138206 |
| Scaffolds_221 | 6137 | 39.4235466536395 |
| Scaffolds_1039 | 6137 | 39.4235466536395 |
| Scaffolds_1027 | 6120 | 56.3847158719791 |
| Scaffolds_492 | 6116 | 39.3919581562602 |
| Scaffolds_52 | 6109 | 39.8167839031572 |
| Scaffolds_866 | 6088 | 56.6316480630335 |
| Scaffolds_612 | 6072 | 35.2369980250165 |
| Scaffolds_798 | 6058 | 36.6545694490267 |
| Scaffolds_270 | 6046 | 65.9504132231405 |
| Scaffolds_604 | 6023 | 39.5719263315082 |
| Scaffolds_946 | 6005 | 30.7871526044267 |
| Scaffolds_862 | 6001 | 60.5328892589509 |
| Scaffolds_429 | 5995 | 47.2745457576263 |
| Scaffolds_872 | 5994 | 60.7702567522508 |
| Scaffolds_295 | 5970 | 39.7221292266488 |
| Scaffolds_561 | 5969 | 59.5680562531391 |
| Scaffolds_401 | 5958 | 35.2566252935257 |
| Scaffolds_618 | 5950 | 50.5542492442056 |
| Scaffolds_266 | 5944 | 45.3093476798924 |
| Scaffolds_804 | 5943 | 43.0973600134522 |
| Scaffolds_631 | 5942 | 58.9976454759502 |
| Scaffolds_756 | 5917 | 38.2198952879581 |
| Scaffolds_522 | 5911 | 49.8055790363483 |
| Scaffolds_465 | 5902 | 64.5106671181849 |
| Scaffolds_649 | 5876 | 43.7755102040816 |
| Scaffolds_376 | 5867 | 46.6019417475728 |
| Scaffolds_1003 | 5862 | 72.4514149335152 |
| Scaffolds_540 | 5858 | 62.0777891504606 |
| Scaffolds_363 | 5842 | 59.3226137529935 |
| Scaffolds_690 | 5833 | 43.0529381531609 |
| Scaffolds_343 | 5826 | 56.4322469982847 |
| Scaffolds_311 | 5813 | 65.7211621110538 |
| Scaffolds_709 | 5803 | 49.5264336146031 |
| Scaffolds_937 | 5783 | 36.9967167789874 |
| Scaffolds_813 | 5780 | 63.0750605326877 |
| Scaffolds_202 | 5771 | 33.8874458874459 |
| Scaffolds_729 | 5771 | 33.8874458874459 |
| Scaffolds_638 | 5766 | 33.3795493934142 |
| Scaffolds_927 | 5742 | 38.722589627567 |
| Scaffolds_498 | 5741 | 71.2793733681462 |
| Scaffolds_530 | 5738 | 37.861372344131 |
| Scaffolds_1028 | 5729 | 54.7880690737834 |
| Scaffolds_289 | 5724 | 46.0719273743017 |
| Scaffolds_1022 | 5704 | 61.8121275849982 |
| Scaffolds_254 | 5694 | 37.048087048087 |
| Scaffolds_678 | 5694 | 37.048087048087 |
| Scaffolds_1026 | 5678 | 58.2365364308342 |
| Scaffolds_774 | 5668 | 41.1847672778561 |
| Scaffolds_464 | 5666 | 36.2433862433862 |
| Scaffolds_654 | 5641 | 39.096545615589 |
| Scaffolds_777 | 5608 | 34.7469707769066 |
| Scaffolds_596 | 5599 | 43.9228984472604 |
| Scaffolds_659 | 5578 | 43.8552490146901 |
| Scaffolds_350 | 5574 | 58.3183936894944 |
| Scaffolds_607 | 5569 | 38.6685806567378 |
| Scaffolds_892 | 5556 | 34.2805755395683 |
| Scaffolds_398 | 5519 | 51.0773130544994 |
| Scaffolds_389 | 5509 | 63.8490839833122 |
| Scaffolds_1041 | 5488 | 70.4297159504734 |
| Scaffolds_688 | 5470 | 39.2583120204604 |
| Scaffolds_528 | 5416 | 62.029520295203 |
| Scaffolds_472 | 5413 | 65.7744138822226 |
| Scaffolds_978 | 5410 | 52.4565940155153 |
| Scaffolds_285 | 5369 | 35.4178298901917 |
| Scaffolds_942 | 5355 | 63.7245754805001 |
| Scaffolds_891 | 5354 | 48.3575961179545 |
| Scaffolds_571 | 5347 | 34.3487198654457 |
| Scaffolds_831 | 5329 | 62.1038814925933 |
| Scaffolds_585 | 5327 | 61.114237478897 |
| Scaffolds_865 | 5319 | 60.229194063498 |
| Scaffolds_837 | 5309 | 36.6647844908714 |
| Scaffolds_466 | 5262 | 33.4219521458412 |
| Scaffolds_906 | 5249 | 34.9324195697697 |
| Scaffolds_725 | 5247 | 37.2690916015997 |
| Scaffolds_441 | 5245 | 69.8037721470756 |
| Scaffolds_967 | 5223 | 62.9473684210526 |
| Scaffolds_129 | 5213 | 39.927161203757 |
| Scaffolds_968 | 5200 | 58.0899308224443 |
| Scaffolds_856 | 5195 | 44.8355452971725 |
| Scaffolds_382 | 5127 | 56.9089846033912 |
| Scaffolds_707 | 5121 | 36.7414634146341 |
| Scaffolds_325 | 5119 | 38.7077884052313 |
| Scaffolds_959 | 5100 | 38.6951410658307 |
| Scaffolds_337 | 5097 | 56.0478337580867 |
| Scaffolds_679 | 5096 | 65.2352941176471 |
| Scaffolds_387 | 5092 | 40.7378335949765 |
| Scaffolds_552 | 5088 | 58.9355852317361 |
| Scaffolds_995 | 5084 | 35.8294025157233 |
| Scaffolds_369 | 5082 | 59.4966574911522 |
| Scaffolds_740 | 5063 | 61.9498717189659 |
| Scaffolds_857 | 5053 | 53.5693098675104 |
| Scaffolds_932 | 5043 | 67.4261937784823 |
| Scaffolds_308 | 5041 | 38.2755203171457 |
| Scaffolds_704 | 5017 | 47.858992232623 |
| Scaffolds_341 | 5010 | 61.268448344635 |
| Scaffolds_620 | 5005 | 34.6176881613096 |
| Scaffolds_1034 | 4989 | 34.7686761466053 |
| Scaffolds_512 | 4977 | 54.5272033728167 |
| Scaffolds_1015 | 4975 | 38.5619602329785 |
| Scaffolds_907 | 4974 | 39.0719164323021 |
| Scaffolds_780 | 4968 | 56.2751407884151 |
| Scaffolds_281 | 4964 | 38.5064412238325 |
| Scaffolds_600 | 4957 | 46.059262245515 |
| Scaffolds_430 | 4946 | 39.6767676767677 |
| Scaffolds_945 | 4938 | 61.3314447592068 |
| Scaffolds_587 | 4929 | 38.7188323535374 |
| Scaffolds_420 | 4928 | 38.0373073803731 |
| Scaffolds_721 | 4921 | 57.6649746192893 |
| Scaffolds_453 | 4916 | 46.869918699187 |
| Scaffolds_715 | 4915 | 65.7653994714373 |
| Scaffolds_821 | 4886 | 64.7239263803681 |
| Scaffolds_474 | 4870 | 61.3254000820681 |
| Scaffolds_263 | 4821 | 58.2590673575129 |
| Scaffolds_1021 | 4819 | 35.7868546547792 |
| Scaffolds_997 | 4797 | 36.6798583628411 |
| Scaffolds_954 | 4789 | 48.7586063008554 |
| Scaffolds_617 | 4777 | 35.8711566617862 |
| Scaffolds_423 | 4775 | 64.1975308641975 |
| Scaffolds_753 | 4768 | 56.8105616093881 |
| Scaffolds_629 | 4743 | 63.2820728881399 |
| Scaffolds_897 | 4742 | 48.7779182469448 |
| Scaffolds_943 | 4742 | 60.9355246523388 |
| Scaffolds_579 | 4735 | 33.5513821481325 |
| Scaffolds_661 | 4730 | 61.8293198141107 |
| Scaffolds_312 | 4718 | 38.56416772554 |
| Scaffolds_783 | 4717 | 49.4174962931582 |
| Scaffolds_958 | 4696 | 55.8085106382979 |
| Scaffolds_439 | 4691 | 35.356762513312 |
| Scaffolds_502 | 4631 | 35.210355987055 |
| Scaffolds_1002 | 4629 | 37.4703216058709 |
| Scaffolds_336 | 4618 | 35.1579402855907 |
| Scaffolds_662 | 4607 | 59.5966167859467 |
| Scaffolds_468 | 4575 | 38.9823105481546 |
| Scaffolds_793 | 4570 | 56.4057717533887 |
| Scaffolds_999 | 4570 | 63.5111499781373 |
| Scaffolds_960 | 4568 | 44.750656167979 |
| Scaffolds_483 | 4567 | 34.9376504047254 |
| Scaffolds_405 | 4556 | 62.8728070175439 |
| Scaffolds_770 | 4520 | 40.3625110521662 |
| Scaffolds_613 | 4513 | 51.162275846801 |
| Scaffolds_816 | 4512 | 33.9902568644818 |
| Scaffolds_938 | 4508 | 38.3200354609929 |
| Scaffolds_795 | 4481 | 59.108138238573 |
| Scaffolds_641 | 4468 | 43.4704830053667 |
| Scaffolds_543 | 4456 | 62.085201793722 |
| Scaffolds_1001 | 4445 | 63.2951224994381 |
| Scaffolds_409 | 4440 | 64.041404140414 |
| Scaffolds_791 | 4425 | 39.3994129600361 |
| Scaffolds_871 | 4422 | 52.7338454586534 |
| Scaffolds_703 | 4383 | 32.8014588557101 |
| Scaffolds_284 | 4326 | 41.824480369515 |
| Scaffolds_431 | 4299 | 25.9121543109459 |
| Scaffolds_878 | 4258 | 33.9981229469733 |
| Scaffolds_639 | 4252 | 66.2828947368421 |
| Scaffolds_806 | 4206 | 63.8242280285036 |
| Scaffolds_473 | 4112 | 44.8007774538387 |
| Scaffolds_432 | 4046 | 65.5555555555556 |
| Scaffolds_859 | 3989 | 61.8081642875031 |
| Scaffolds_294 | 3624 | 75 |
| Scaffolds_815 | 3444 | 33.4106728538283 |
| Scaffolds_970 | 3394 | 53.0017657445556 |
| Scaffolds_818 | 3371 | 65.8167803142603 |
| Scaffolds_42 | 3169 | 35.4238890639773 |
| Scaffolds_838 | 3122 | 23.320537428023 |
| Scaffolds_841 | 2616 | 26.2996941896024 |
