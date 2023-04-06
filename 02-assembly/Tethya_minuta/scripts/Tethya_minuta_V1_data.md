## *Tethya minuta* GW41610

## Nanopore reads

 Feb 11 18:21 Tminuta_MinION_nanopore_1D_combined.fastq

## Illumina sequences
*  Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq
*  Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq
## RNA
* 1link_Trimmomatic_paired_GW41610_TGCACGAT_1.fastq
* 1link_Trimmomatic_paired_GW41610_TGCACGAT_2.fastq
## WTDBG2

```shell
#!/bin/bash
#
#SBATCH --job-name=WTDBG
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000
 
wtdbg2 -x ont -t 112 -L 500  -i ../Tminuta_MinION_nanopore_1D_combined.fastq -fo Tminuta_assembly_WTDBG_nanopore

wtpoa-cns -t 112 -i Tminuta_assembly_WTDBG_nanopore.ctg.lay.gz -fo Tminuta_assembly_WTDBG_nanopore.ctg.fa 

minimap2 -t 112 -x map-pb -a Tminuta_assembly_WTDBG_nanopore.ctg.fa \
../Tminuta_MinION_nanopore_1D_combined.fastq | samtools view -Sb - >Tminuta_assembly_WTDBG_nanopore.ctg.map.bam

samtools sort Tminuta_assembly_WTDBG_nanopore.ctg.map.bam  >Tminuta_assembly_WTDBG_nanopore.ctg.map.srt.bam

samtools view Tminuta_assembly_WTDBG_nanopore.ctg.map.srt.bam | wtpoa-cns -t 112 -d Tminuta_assembly_WTDBG_nanopore.ctg.fa -i - \
-fo Tminuta_assembly_WTDBG_nanopore.ctg.2nd.fa

bwa index Tminuta_assembly_WTDBG_nanopore.ctg.2nd.fa

bwa mem -t 112 Tminuta_assembly_WTDBG_nanopore.ctg.2nd.fa \
Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq | samtools view -Sb - >sr.bam
samtools sort sr.bam >sr.srt.bam
samtools view sr.srt.bam | wtpoa-cns -t 112 -x sam-sr -d Tminuta_assembly_WTDBG_nanopore.ctg.fa -i - -fo Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa
```


# statistics
```shell
perl sequence_counter.pl Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa 

longest scaffold = 3190713
shortest scaffold = 3434
N50 = 133218
number of seuqneces = 2193
number of N = 0
N percent = 0
total length = 126023507

```
# Busco 
```shell
# BUSCO version is: 3.1.0 
# The lineage dataset is: metazoa_odb9 (Creation date: 2016-02-13, number of species: 65, number of BUSCOs: 978)
# To reproduce this run: python ../scripts/run_BUSCO.py -i Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa -o OUTPUT -l /naslx/projects/pn69xe/di52zuc/program/busco-master/metazoa_odb9/ -m genome -c 28 -sp Tethya_wilh
#
# Summarized benchmarking in BUSCO notation for file Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa
# BUSCO was run in mode: genome

        C:66.1%[S:62.9%,D:3.2%],F:9.3%,M:24.6%,n:978

        646     Complete BUSCOs (C)
        615     Complete and single-copy BUSCOs (S)
        31      Complete and duplicated BUSCOs (D)
        91      Fragmented BUSCOs (F)
        241     Missing BUSCOs (M)
        978     Total BUSCO groups searched

```


## Read mapping
```shell
690421 + 0 in total (QC-passed reads + QC-failed reads)
77365 + 0 secondary
67564 + 0 supplementary
0 + 0 duplicates
523629 + 0 mapped (75.84% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

# This file was produced by samtools stats (1.9+htslib-1.9) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats sr.bam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     2c63fcf1        18a38391        237dc3af
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    42688960
SN      filtered sequences:     0
SN      sequences:      42688960
SN      is sorted:      0
SN      1st fragments:  21344480
SN      last fragments: 21344480
SN      reads mapped:   35014289
SN      reads mapped and paired:        33620028        # paired-end technology bit set + both mates mapped
SN      reads unmapped: 7674671
SN      reads properly paired:  33038580        # proper-pair bit set
SN      reads paired:   42688960        # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      748337  # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 0
SN      total length:   3071560344      # ignores clipping
SN      total first fragment length:    1673113715      # ignores clipping
SN      total last fragment length:     1398446629      # ignores clipping
SN      bases mapped:   2614723274      # ignores clipping
SN      bases mapped (cigar):   2476889872      # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     74683482        # from NM fields
SN      error rate:     3.015212e-02    # mismatches / bases mapped (cigar)
SN      average length: 71
SN      average first fragment length:  78
SN      average last fragment length:   66
SN      maximum length: 93
SN      maximum first fragment length:  93
SN      maximum last fragment length:   93
SN      average quality:        36.9
SN      insert size average:    163.5
SN      insert size standard deviation: 78.4
SN      inward oriented pairs:  14018154
SN      outward oriented pairs: 2558911
SN      pairs with other orientation:   11065
SN      pairs on different chromosomes: 229916
SN      percentage of properly paired reads (%):        77.4
```

```shell
#!/bin/bash
#
#SBATCH --job-name=Tophat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000

building Bowtie2 index for Tethya minuta nanopore canu assembly

bowtie2-build --large-index  minuta.contigs.fa minuta.contig 

mapping using Tophat2

tophat2 --no-mixed minuta.contig link_Trimmomatic_paired_GW41610_TGCACGAT_1.fastq,link_Trimmomatic_paired_GW41610_TGCACGAT_2.fastq


```

```shell
# This file was produced by samtools stats (1.9+htslib-1.9) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats accepted_hits.bam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     4973df2a        07e1816a        42186265
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    11045752
SN      filtered sequences:     0
SN      sequences:      11045752
SN      is sorted:      1
SN      1st fragments:  11045752
SN      last fragments: 0
SN      reads mapped:   11045752
SN      reads mapped and paired:        0       # paired-end technology bit set + both mates mapped
SN      reads unmapped: 0
SN      reads properly paired:  0       # proper-pair bit set
SN      reads paired:   0       # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      46718   # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 1028441
SN      total length:   1014334906      # ignores clipping
SN      total first fragment length:    1014334906      # ignores clipping
SN      total last fragment length:     0       # ignores clipping
SN      bases mapped:   1014334906      # ignores clipping
SN      bases mapped (cigar):   1014334906      # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     9790556 # from NM fields
SN      error rate:     9.652193e-03    # mismatches / bases mapped (cigar)
SN      average length: 91
SN      average first fragment length:  92
SN      average last fragment length:   0
SN      maximum length: 92
SN      maximum first fragment length:  0
SN      maximum last fragment length:   0
SN      average quality:        37.5
SN      insert size average:    0.0
```
## SSPACE_Longread Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa
Nanopore reads: 4.2G Tminuta_MinION_nanopore_1D_combined.fastq

```shell
perl sequence_counter.pl T-minuta_nanopore_scaffolded.fasta
longest scaffold = 3190713
shortest scaffold = 3434
N50 = 183430
number of sequences = 1690
number of N|n-position = 1315559
N|n percent = 1.0331150065134
total length = 127339066

```
## RNA mapping Hisat2

```shell
#!/bin/bash
#
#SBATCH --job-name=Hisat2
#SBATCH --ntasks=32
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=Hisat2.err
#SBATCH --mem=125G

# build index for query sequence for Hisat2

hisat2-build T-minuta_nanopore_scaffolded.fasta \
scaffolds.fasta 

# run HISAT v. 2.1.0 

hisat2 -q --no-mixed --no-unal --no-spliced-alignment --threads 32 -x scaffolds.fasta \
-1 1link_Trimmomatic_paired_GW41610_TGCACGAT_1.fastq \
-2 2link_Trimmomatic_paired_GW41610_TGCACGAT_2.fastq \
-S scaffolds.sam 2>errorlog.txt

#writing file containing stats

samtools stats scaffolds.sam >bamstats.txt
samtools flagstat scaffolds.sam >flagstats.txt

```

# samtools flagstat

```shell
13755628 + 0 in total (QC-passed reads + QC-failed reads)
370782 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
13755628 + 0 mapped (100.00% : N/A)
13384846 + 0 paired in sequencing
6692423 + 0 read1
6692423 + 0 read2
12734152 + 0 properly paired (95.14% : N/A)
13384846 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
57324 + 0 with mate mapped to a different chr
57324 + 0 with mate mapped to a different chr (mapQ>=5)
```
# samtools stats
```shell
# This file was produced by samtools stats (1.8+htslib-1.8) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats scaffolds.sam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     b7fc365a        cc0a1057        c8dcccd3
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    13384846
SN      filtered sequences:     0
SN      sequences:      13384846
SN      is sorted:      0
SN      1st fragments:  6692423
SN      last fragments: 6692423
SN      reads mapped:   13384846
SN      reads mapped and paired:        13384846        # paired-end technology bit set + both mates mapped
SN      reads unmapped: 0
SN      reads properly paired:  12734152        # proper-pair bit set
SN      reads paired:   13384846        # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      131896  # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 370782
SN      total length:   1229095327      # ignores clipping
SN      bases mapped:   1229095327      # ignores clipping
SN      bases mapped (cigar):   1208095084      # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     4208612 # from NM fields
SN      error rate:     3.483676e-03    # mismatches / bases mapped (cigar)
SN      average length: 91
SN      maximum length: 92
SN      average quality:        37.5
SN      insert size average:    233.0
SN      insert size standard deviation: 312.4
SN      inward oriented pairs:  6276710
SN      outward oriented pairs: 294943
SN      pairs with other orientation:   1730
SN      pairs on different chromosomes: 28662

```
## GapCloser
```
#!/bin/bash
#
#SBATCH --job-name=GapCloser
#SBATCH --ntasks=64
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --err=GapCloser_wilhelma.err
#SBATCH --mem=300G


GapCloser -a T-minuta_scaffolded_nanopore.fasta \
-b Library_GapCloser.txt -o gap_closed_T-minuta_scaffolded_nanopore.fasta -t 64
```
```bash
#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=200
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=50
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_1.fastq
q2=Trimmomatic_paired_L-025_FD_AGAGTAGA-TAAGGCGA_2.fastq

```
