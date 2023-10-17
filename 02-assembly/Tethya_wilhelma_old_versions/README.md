# *T. wilhelma* assemblies #

| Version | # scaffolds | Total size (Mb) | N50 (Mb) | N gaps (kb) |
| --- | --: | --: | :-- | --: |
| v1 | 7947 | 145.36 | 0.065 | 1539.6 | 
| v2 | 1353 | 139.49 | 5.5 | 2252.2 | 
| v3 | 967 | 138.92 | 6.1 | 1069.9 | 
| v4 | 891 | 138.93 | 6.7 | 1077.5 |

```
supernova run --id=Twilhelma_final --fastqs=~/Tethya_wilhelma_10x/all_data --sample=TW --maxreads="all" --localcores=80 --localmem=755 --accept-extreme-coverage


hicstuff iteralign -g gaps_closed_scaffolded_10_x_raw.fasta  -t 28 -o for.bam output_DTG-HiC-396_R1_001_forward_paired.fq.gz

samtools sort -@ 28 -O bam -T for.sort -n -o for.sorted.bam for.bam

hicstuff iteralign -g gaps_closed_scaffolded_10_x_raw.fasta  -t 28 -o rev.bam output_DTG-HiC-396_R2_001_reverse_paired.fq.gz

samtools sort -@ 28 -O bam -T rev.sort -n -o rev.sorted.bam rev.bam

hicstuff pipeline -e DpnII -g gaps_closed_scaffolded_10_x_raw.fasta -S bam -n -o .  -t 28 for.sorted.bam rev.sorted.bam


GapCloser -a Dovetail_assembly__sponge_tethya_wilhelma_26Apr2018_OM745.fasta -b Library_GapCloser.txt -o gaps_closed_scaffolded_short_reads.final.scaffolds.fasta -t 40

blastp -db SwissProt -query protein_predicted_genes.gff.fasta -out tabular_blast_wilhelma_swissprot_out.txt -outfmt 6
```

### mapping RNAseq ###

```
hisat2-build coding_regions_predicted_genes_coding_sequences.gff.fasta scaffolds.fasta 
hisat2 -q --no-mixed --no-unal --no-spliced-alignment --threads 64 -x scaffolds.fasta -1 T-wilhelma_Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq -2 T-wilhelma_Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq -S scaffolds.sam 2> errorlog.txt
samtools flagstat -@64  scaffolds.sam > flagstats.txt

hisat2 -q --no-mixed --no-unal --no-spliced-alignment --threads 64 -x scaffolds.fasta -1 T-wilhelma_Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq -2 T-wilhelma_Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq -S scaffolds.sam 2> errorlog.txt
samtools stats -@64  scaffolds.sam > samstats.txt
```

### flagstats.txt

```shell
91726772 + 0 in total (QC-passed reads + QC-failed reads)
2854186 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
91726772 + 0 mapped (100.00% : N/A)
88872586 + 0 paired in sequencing
44436293 + 0 read1
44436293 + 0 read2
87613378 + 0 properly paired (98.58% : N/A)
88872586 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
614614 + 0 with mate mapped to a different chr
614614 + 0 with mate mapped to a different chr (mapQ>=5)

```

### samstats.txt
```shell
#samtstats
# This file was produced by samtools stats (1.8+htslib-1.8) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats -@64 scaffolds.sam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     2c25e20c        e58d8c69        9131eda9
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    88872586
SN      filtered sequences:     0
SN      sequences:      88872586
SN      is sorted:      0
SN      1st fragments:  44436293
SN      last fragments: 44436293
SN      reads mapped:   88872586
SN      reads mapped and paired:        88872586        # paired-end technology bit set + both mates mapped
SN      reads unmapped: 0
SN      reads properly paired:  87613378        # proper-pair bit set
SN      reads paired:   88872586        # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      1251492 # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 2854186
SN      total length:   10398092562     # ignores clipping
SN      bases mapped:   10398092562     # ignores clipping
SN      bases mapped (cigar):   10352941249     # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     46182328        # from NM fields
SN      error rate:     4.460793e-03    # mismatches / bases mapped (cigar)
SN      average length: 117
SN      maximum length: 117
SN      average quality:        35.5
SN      insert size average:    194.9
SN      insert size standard deviation: 39.0
SN      inward oriented pairs:  43750967
SN      outward oriented pairs: 322038
SN      pairs with other orientation:   259
SN      pairs on different chromosomes: 307307
```

