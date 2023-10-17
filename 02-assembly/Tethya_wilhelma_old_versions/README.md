# *T. wilhelma* assemblies #

| Version | # scaffolds | Total size (Mb) | N50 (Mb) | N gaps (kb) |
| --- | --: | --: | :-- | --: |
| v1 | 7947 | 145.36 | 0.065 | 1539.6 | 
| v2 | 1353 | 139.49 | 5.5 | 2252.2 | 
| v3 | 967 | 138.92 | 6.1 | 1069.9 | 
| v4 | 891 | 138.93 | 6.7 | 1077.5 |

# v4 #

```
~/git/genomeGTFtools/repeat2gtf.py -l Tethya_wilhelma_V4_P_RNA_scaffold.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.n_gaps.gff
# download SRA data and map
# to get coverage and identify bacterial contigs

# genomic DNA paired end reads
# https://www.ncbi.nlm.nih.gov/sra/SRX1149776[accn]
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --split-files --gzip SRR2163223

# genomic mate pair
# https://www.ncbi.nlm.nih.gov/sra/SRX1149801[accn]
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --split-files --gzip SRR2296844

# RNAseq paired end reads strand specific library prep
# https://www.ncbi.nlm.nih.gov/sra/SRX2175444[accn]
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --split-files --gzip SRR4255675

# metatranscriptomic paired reads
# https://www.ncbi.nlm.nih.gov/sra/ERX9587876[accn]
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --split-files --gzip ERR10048047
#Read 28769610 spots for ERR10048047
#Written 28769610 spots for ERR10048047


# convert bw and bed with
# https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

# genomic moleculo long reads
# https://www.ncbi.nlm.nih.gov/sra/SRX2665261[accn]
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip SRR5369934
#Read 125150 spots for SRR5369934
#Written 125150 spots for SRR5369934

~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --split-files --gzip ERR10048047
{ ~/hisat2-2.2.1/hisat2 -q -x Tethya_wilhelma_V4_P_RNA_scaffold.fasta -1 ERR10048047_1.fastq.gz -2 ERR10048047_2.fastq.gz -p 6 --dta --max-intronlen 20000 2>&3 | ~/samtools-1.14/samtools sort - -o ERR10048047_vs_V4_hisat2.bam; } 3> ERR10048047_vs_V4_hisat2.log
~/samtools-1.14/samtools view ERR10048047_vs_V4_hisat2.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > ERR10048047_vs_V4_hisat2.hits_from_bam.txt
~/git/lavaLampPlot/hits_to_coverage.py -b ERR10048047_vs_V4_hisat2.hits_from_bam.txt -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta -l 75 > ERR10048047_vs_V4_hisat2.gc_cov.tab

# map to genome
~/minimap2-2.23_x64-linux/minimap2 -a -x map-hifi --secondary=no Tethya_wilhelma_V4_P_RNA_scaffold.fasta SRR5369934.fastq.gz | ~/samtools-1.14/samtools sort - -o Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam
~/samtools-1.14/samtools index Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam
~/samtools-1.14/samtools flagstat Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam
#146165 + 0 in total (QC-passed reads + QC-failed reads)
#125150 + 0 primary
#0 + 0 secondary
#21015 + 0 supplementary
#139704 + 0 mapped (95.58% : N/A)
#118689 + 0 primary mapped (94.84% : N/A)
bedtools genomecov -ibam Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam -bg > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bed

~/git/lavaLampPlot/hits_to_coverage.py -g Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bed -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta -l 1 > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.gc_cov.tab
# IGNORING scaffold721: len=2563 GC=40.19 COVERAGE=0.00
# IGNORING scaffold727: len=2497 GC=31.56 COVERAGE=0.00
# IGNORING scaffold801: len=2010 GC=40.75 COVERAGE=0.00
# IGNORING scaffold819: len=1915 GC=57.18 COVERAGE=0.00
# IGNORING scaffold820: len=1914 GC=39.34 COVERAGE=0.00
# IGNORING scaffold838: len=1808 GC=38.05 COVERAGE=0.00
# IGNORING scaffold839: len=1805 GC=40.66 COVERAGE=0.00
# IGNORING scaffold888: len=1603 GC=65.07 COVERAGE=0.00
# IGNORING scaffold921: len=1402 GC=36.45 COVERAGE=0.00
# IGNORING scaffold853: len=1770 GC=40.11 COVERAGE=0.00
# IGNORING scaffold937: len=1197 GC=57.56 COVERAGE=0.00
# IGNORING scaffold939: len=1187 GC=38.42 COVERAGE=0.00
# IGNORING scaffold941: len=1161 GC=41.69 COVERAGE=0.00
# IGNORING scaffold952: len=1068 GC=53.37 COVERAGE=0.00
# IGNORING scaffold953: len=1062 GC=35.97 COVERAGE=0.00
# IGNORING scaffold954: len=1058 GC=39.79 COVERAGE=0.00
# IGNORING scaffold959: len=1016 GC=36.12 COVERAGE=0.00
# IGNORING scaffold960: len=1007 GC=38.23 COVERAGE=0.00
```

# v4 binning bacteria #

```
~/samtools-1.14/samtools view Anoxia_shock_combined_v4.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > Anoxia_shock_combined_v4.hits_from_bam.txt

~/git/lavaLampPlot/hits_to_coverage.py -g Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bed -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta --rna-read-length 50 -r Anoxia_shock_combined_v4.hits_from_bam.txt > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.gc_cov.tab

~/git/lavaLampPlot/hits_to_coverage.py -g Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.bed.gz -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta --rna-read-length 50 -r Anoxia_shock_combined_v4.hits_from_bam.txt > Tethya_wilhelma_V4_P_RNA_scaffold.hisat2_dna.gc_cov.tab

cut -f 1 Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.tab > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names
getAinB.py Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names Tethya_wilhelma_V4_P_RNA_scaffold.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.b47.fasta

```

# v4 genome browser including bacterial contigs #

```
~/samtools-1.14/samtools faidx Tethya_wilhelma_V4_P_RNA_scaffold.fasta

ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.fasta
ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.fasta.fai 
ln -s /user/Tethya_wilhelma_V4/predicted_genes_Tethya_wilhelma_V4_below_50_perc.no_comment.gff
/var/www/html/jbrowse/bin/prepare-refseqs.pl --indexed_fasta Tethya_wilhelma_V4_P_RNA_scaffold.fasta --out ./
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff predicted_genes_Tethya_wilhelma_V4_below_50_perc.no_comment.gff --trackType CanvasFeatures --trackLabel AUGUSTUS --out ./

ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam
ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bam.bai 

ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.twiV1_curated_minimap2.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff Tethya_wilhelma_V4_P_RNA_scaffold.twiV1_curated_minimap2.gff --trackType CanvasFeatures --trackLabel v1-filtered --out ./
ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.twiV1_stringtie_minimap2.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff Tethya_wilhelma_V4_P_RNA_scaffold.twiV1_stringtie_minimap2.gff --trackType CanvasFeatures --trackLabel v1-stringTie --out ./

ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.twilhelma_trinity_ss_norm_minimap2.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff Tethya_wilhelma_V4_P_RNA_scaffold.twilhelma_trinity_ss_norm_minimap2.gff --trackType CanvasFeatures --trackLabel Trinity-SS-minimap --out ./

ln -s /user/Tethya_wilhelma_V4/predicted_genes_Tethya_wilhelma_V4_below_50_perc.prot.vs_aque.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff predicted_genes_Tethya_wilhelma_V4_below_50_perc.prot.vs_aque.gff --trackType CanvasFeatures --trackLabel blast-vs-aque --out ./
ln -s /user/Tethya_wilhelma_V4/predicted_genes_Tethya_wilhelma_V4_below_50_perc.prot.vs_model.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff predicted_genes_Tethya_wilhelma_V4_below_50_perc.prot.vs_model.gff --trackType CanvasFeatures --trackLabel blast-vs-models --out ./

ln -s /user/Tethya_wilhelma_V4/Tethya_wilhelma_V4_P_RNA_scaffold.n_gaps.gff
/var/www/html/jbrowse/bin/flatfile-to-json.pl --gff Tethya_wilhelma_V4_P_RNA_scaffold.n_gaps.gff --trackType CanvasFeatures --trackLabel N-gaps --out ./

```

# Twi v3 #
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

