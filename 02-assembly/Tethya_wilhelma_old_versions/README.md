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

sizecutter.py -f Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes
~/git/genomeGTFtools/rename_gtf_contigs.py -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -n -g Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.bed.gz > Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed
# Reading conversion file Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector
# Found names for 557 contigs
# Counted 98593420 lines
# Converted 91299383 lines and could not change 7294037
~/ucsc_genome_tools/bedSort Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed
~/ucsc_genome_tools/bedGraphToBigWig Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bw

~/git/genomeGTFtools/rename_gtf_contigs.py -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -n -g Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.bed.gz > Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed
# Reading conversion file Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector
# Found names for 557 contigs
# Counted 27471719 lines
# Converted 27443244 lines and could not change 28475
~/ucsc_genome_tools/bedSort Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed
~/ucsc_genome_tools/bedGraphToBigWig Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes Twi_RNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.rnum.bw


```

### RNAseq mapping 

```shell
237,305,498 + 0 in total (QC-passed reads + QC-failed reads)
39,365,382 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
224,608,585 + 0 mapped (94.65% : N/A)
197,940,116 + 0 paired in sequencing
98,970,058 + 0 read1
98,970,058 + 0 read2
169,714,632 + 0 properly paired (85.74% : N/A)
183,839,570 + 0 with itself and mate mapped
1403633 + 0 singletons (0.71% : N/A)
1921096 + 0 with mate mapped to a different chr
1011724 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Busco scores
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

### v4 binning bacteria ###

```
~/samtools-1.14/samtools view Anoxia_shock_combined_v4.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > Anoxia_shock_combined_v4.hits_from_bam.txt

~/git/lavaLampPlot/hits_to_coverage.py -g Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.bed -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta --rna-read-length 50 -r Anoxia_shock_combined_v4.hits_from_bam.txt > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.gc_cov.tab

~/git/lavaLampPlot/hits_to_coverage.py -g Twi_DNA_vs_Tethya_wilhelma_V4_P_RNA_scaffold.bed.gz -f Tethya_wilhelma_V4_P_RNA_scaffold.fasta --rna-read-length 50 -r Anoxia_shock_combined_v4.hits_from_bam.txt > Tethya_wilhelma_V4_P_RNA_scaffold.hisat2_dna.gc_cov.tab

cut -f 1 Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.tab > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names
getAinB.py Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names Tethya_wilhelma_V4_P_RNA_scaffold.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.b47.fasta

```

### v4 genome browser including bacterial contigs ###

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

### v2 to v3 scaffolding

```shell
#!/bin/bash

export PATH=$PATH:/user/HaploMerger2_20180603/lastz_1.02.00_unbuntu64bit/ 
export PATH=$PATH:/user/HaploMerger2_20180603/chainNet_jksrc20100603_ubuntu64bit/  
export PATH=$PATH:/user/exonerate-2.2.0-x86_64/bin/
#### ===========================================================
#### A1. Initiate data to be ready for use in the next steps.
#### ===========================================================
rm -f -r Twilhelma.seq Twilhelmax.seq Twilhelma.sizes Twilhelmax.sizes Twilhelmax.fa.gz 
perl ../bin/initiation.pl --faSplit --faToNib --faSize --Species Twilhelma --Force --Delete \
        1>_A1.initiation.log 2>>_A1.initiation.log
ln -s Twilhelma.fa.gz Twilhelmax.fa.gz
ln -s Twilhelma.sizes Twilhelmax.sizes
ln -s Twilhelma.seq Twilhelmax.seq
#### ===========================================================
#### run all-against-all whole genome alignment
#### ===========================================================
rm -f -r Twilhelma.Twilhelmax.result/raw.axt
perl ../bin/HM_all_lastz_mThreads.pl --Species Twilhelma Twilhelmax --noself --threads=112 --identity=80 \
  --targetSize=13000000 --querySize=1300000000 --Force --Delete \
  1>_A1.all_lastz.log 2>>_A1.all_lastz.log
#### ===========================================================  
#### A2 compute reciprocally-best alignment
#### ===========================================================
perl ../bin/HM_axtChainRecipBestNet.pl --rbestNet --axtChain --tbest --zeroMinSpaceNet \
--threads=112 \
--axtSuffix=axt.noself --linearGap=medium --minScore=5000 --minSpace=100 --minScore2=10000 \
--Species Twilhelma Twilhelmax --Force --Delete \
1>_A2.axtChainRecipBestNet.log 2>>_A2.axtChainRecipBestNet.log
#### ===========================================================
#### A3
#### ===========================================================
perl ../bin/HM_pathFinder_preparation.pl --Species Twilhelma Twilhelmax \
  --scoreScheme=$scoreScheme --filter=100000 --Force --Delete \
  1>_A3.pathFinder_preparation.log 2>>_A3.pathFinder_preparation.log
perl ../bin/HM_pathFinder.pl --Species Twilhelma Twilhelmax \
--Force --Delete --scoreScheme=score --filter=200000 \
--NsLCsFilter=90 --noSelfLoop=1 --noStrandConflict=1 --breakingMode=2 \
--misjoin_aliFilter=5000000 --misjoin_overhangFilter=50000 --escapeFilter=100 \
1>_A3.pathFinder.log 2>>_A3.pathFinder.log
perl ../bin/XHM_remove_misjoin_scf.pl --scf=Twilhelma.fa.gz --misjoin=Twilhelma.Twilhelmax.result/hm.assembly_errors \
--aliFilter=5000000 --overhangFilter=50000 2>_A3.misjoin_processing.log \
| perl ../bin/faDnaPolishing.pl --legalizing --noLeadingN --removeShortSeq=200 2>_A3.faDnaPolishing.log \
| gzip -c >Twilhelma_A.fa.gz
#### ===========================================================  
#### B1
#### ===========================================================
rm -f -r Twilhelma_A.seq Twilhelma_Ax.seq Twilhelma_A.sizes Twilhelma_Ax.sizes Twilhelma_Ax.fa.gz
perl ../bin/initiation.pl --faSplit --faToNib --faSize --Species Twilhelma_A --Force --Delete \
  1>_B1.initiation.log 2>>_B1.initiation.log
ln -s Twilhelma_A.fa.gz Twilhelma_Ax.fa.gz
ln -s Twilhelma_A.sizes Twilhelma_Ax.sizes
ln -s Twilhelma_A.seq Twilhelma_Ax.seq
## require two control files: 1) all_lastz.ctl|$name.{name}x.ctl, and 2) scoreMatrix.q|$name.${name}x.q
rm -f -r Twilhelma_A.Twilhelma_Ax.result/raw.axt
perl ../bin/HM_all_lastz_mThreads.pl --Species Twilhelma_A Twilhelma_Ax --notrivial --radius=5000 --threads=112 --identity=80 \
  --targetSize=13000000 --querySize=1300000000 --Force --Delete \
  1>_B1.all_lastz.log 2>>_B1.all_lastz.log
#### ===========================================================
#### B2 and B3
#### ===========================================================
#### run axtChainRecipBestNet
#### ===========================================================
perl ../bin/HM_axtChainRecipBestNet.pl --rbestNet --axtChain --tbest --zeroMinSpaceNet --threads=112 \
--axtSuffix=axt.self.notrivial --linearGap=medium --minScore=1000 --minScore2=2000 \
--Species Twilhelma_A Twilhelma_Ax --Force --Delete \
1>_B2.axtChainRecipBestNet.log 2>>_B2.axtChainRecipBestNet.log
#### ===========================================================
#### run netToMaf
#### ===========================================================
rm -f Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.gz
ln -s all.rbest.net.gz Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.gz
perl ../bin/HM_netToAxtMaf.pl    \
  --Species Twilhelma_A Twilhelma_Ax       \
  --Force                        \
  --netFile=mafFiltered.net.gz   \
  --chainFile=all.rbest.chain.gz \
  1>_B2.netToMaf.log 2>>_B2.netToMaf.log
perl ../bin/mafSplitByTarget.pl                \
  Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.maf.gz \
  1>>_B2.netToMaf.log 2>>_B2.netToMaf.log
tar -czf Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.maf.tar.gz Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.maf
rm -f -r Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.maf
rm -f Twilhelma_A.Twilhelma_Ax.result/mafFiltered.net.maf.gz
#### ===========================================================
#### preparation - to convert the net infomation
#### ===========================================================
perl ../bin/HM_pathFinder_preparation.pl --Species Twilhelma_A Twilhelma_Ax --Force --Delete \
  --scoreScheme=score --filter=20000 \
  1>_B3.pathFinder_preparation.log 2>>_B3.pathFinder_preparation.log
#### ===========================================================
#### find the paths and detect misjoins
#### ===========================================================
perl ../bin/HM_pathFinder.pl --Species Twilhelma_A Twilhelma_Ax --Force --Delete --scoreScheme=score --filter=200000 \
        --NsLCsFilter=90 --noSelfLoop=1 --noStrandConflict=1 --breakingMode=2 \
        --misjoin_aliFilter=5000000 --misjoin_overhangFilter=50000 --escapeFilter=100 \
        1>_B3.pathFinder.log 2>>_B3.pathFinder.log
#### ===========================================================
#### output haplomic assemblies
#### ===========================================================
perl ../bin/XHM_haploMerger.pl --Species Twilhelma_A Twilhelma_Ax --Force --Delete  \
  --selectLongHaplotype --minOverlap=0  \
  1>_B3.haploMerger.log 2>>_B3.haploMerger.log
#### ===========================================================
#### remove leading and trailing Ns from the scaffold sequences
#### ===========================================================
gunzip -c Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds.fa.gz \
        | perl ../bin/faDnaPolishing.pl --legalizing --noLeadingN --removeShortSeq=500 2>_B3.faDnaPolishing.log \
        | gzip -c > Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_new.fa.gz
gunzip -c Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_alt.fa.gz \
        | perl ../bin/faDnaPolishing.pl --legalizing --noLeadingN --removeShortSeq=500 2>>_B3.faDnaPolishing.log \
        | gzip -c > Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_alt_new.fa.gz
gunzip -c Twilhelma_A.Twilhelma_Ax.result/unpaired.fa.gz \
        | perl ../bin/faDnaPolishing.pl --legalizing --noLeadingN --removeShortSeq=500 2>>_B3.faDnaPolishing.log \
        | gzip -c > Twilhelma_A.Twilhelma_Ax.result/unpaired_new.fa.gz

rm -f Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds.fa.gz \
  Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_alt.fa.gz \
  Twilhelma_A.Twilhelma_Ax.result/unpaired.fa.gz

mv Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_new.fa.gz Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds.fa.gz
mv Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_alt_new.fa.gz Twilhelma_A.Twilhelma_Ax.result/optiNewScaffolds_alt.fa.gz
mv Twilhelma_A.Twilhelma_Ax.result/unpaired_new.fa.gz Twilhelma_A.Twilhelma_Ax.result/unpaired.fa.gz
#### ===========================================================
#### B4
#### ===========================================================
#### refine the unpaired sequences
#### ===========================================================
perl ../bin/HM_unincorpRefiner.pl --Species Twilhelma_A Twilhelma_Ax \
--runLastzChainNet=1 --threads=112 \
--identity=75 --maskFilter=85 \
--redundantFilter=85 --Force --Delete \
1>_B4.unincorpRefiner.log 2>>_B4.unincorpRefiner.log

mv -f _hm.un_initiation.log _B4.un_initiation.log
mv -f _hm.un_all_lastz.log _B4.un_all_lastz.log
mv -f _hm.un_axtChainRecipBestNet.log _B4.un_axtChainRecipBestNet.log

```
```shell
longest scaffold = 18537106
shortest scaffold = 708
N50 = 5547405
Number of sequences = 1291
Number of N-positions= 2247036
total length = 138669404
```

## v3 SSPACE Long read
Using [SSPACE_LongRead](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-211)

```
perl /user/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl \
-c Twilhelma_A_ref.fasta \
-b Twilhelma_moleculo_reads_sspace \
-i 80 -t 80 -o 5000 -l 3 -k 1 \
-p Galaxy173_moleculo_read_T_wilhelma.fastq 2> moleculo_reads_errorlog.txt >moleculo_reads_contigs_logfile.txt
```


```shell
scaffolds.fasta > 
longest scaffold = 27149804
shortest scaffold = 708
N50 = 6150662
Number of sequences = 1035
Number of N-positions = 2391579
total length = 138812862
```

### v3 RNA mapping

```shell

# build index for query sequence for Hisat2

hisat2-build gaps_T-wilhelma_after_HM_and_SSPACE.scaffolds.fasta  Twilhelma


# run HISAT v. 2.1.0 

hisat2 -q --no-mixed --no-unal --no-spliced-alignment --threads 56 \
-x Twilhelma -1 ../Trimmomatic_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq \
-2../Trimmomatic_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq \
-S T_wilhelma_mapped_RNA.sam 2>errorlog.txt

# sort sam to bam

samtools view -F 4 -T Twilhelma.fa -Sb Twilhema_mapped_RNA.sam | samtools sort -@ 56 > Twilhelma_mapped_RNA.bam

```


# *Tethya wilhelma* V2

```shell
longest scaffold = 18535996
shortest scaffold = 410
N50 = 5547405
number of sequences = 1363
number of N|n-position = 2255791
N|n percent = 1.61712446294181
total length = 139493963
```

# v1 #

```shell
longest scaffold = 659.656 
shortest scaffold = 1.004 
N50 = 73.701 
number of sequences = 6.420 
number of N = 1.348.047 
N percent = 1.07268270022062 
total length = 125.670.620
```

```
graphmap align -t 28 -r twilhelma_scaffolds_v1.fasta  -d moleculo_reads.fastq -o mapped_moleculo >logfile 2> errorlog  

720151 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
664445 + 0 mapped (92.26% : N/A)
```

