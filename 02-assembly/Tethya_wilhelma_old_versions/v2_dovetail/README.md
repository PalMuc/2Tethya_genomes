# tethya-v2

```
#!/bin/bash

cd ~/genomes/tethya_wilhelma-genome/rawdata

### download DNA short reads

~/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --split-files --gzip SRR2163223

mv SRR2163223_1.fastq.gz twilhelma_2014_dna_SRR2163223_1.fastq.gz
mv SRR2163223_2.fastq.gz twilhelma_2014_dna_SRR2163223_2.fastq.gz

### download RNAseq reads

~/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --split-files --gzip --defline-seq '@$sn[_$rn]/$ri' SRR4255675

mv SRR4255675_1.fastq.gz twilhelma_rnaseq_SRR4255675_1.fastq.gz
mv SRR4255675_2.fastq.gz twilhelma_rnaseq_SRR4255675_2.fastq.gz


### Tethya V1 DNA MAPPING FOR COVERAGE

hisat2-build twilhelma_scaffolds_v1.fasta twilhelma_scaffolds_v1.fasta

hisat2 -q -x twilhelma_scaffolds_v1.fasta -1 ../rawdata/twilhelma_2014_dna_SRR2163223_1.fastq.gz -2 ../rawdata/twilhelma_2014_dna_SRR2163223_2.fastq.gz -S twilhelma_2014_vs_scaffolds_v1.sam -p 4 --no-spliced-alignment 2> twilhelma_2014_vs_scaffolds_v1.log

~/samtools-1.8/samtools view -Su twilhelma_2014_vs_scaffolds_v1.sam | ~/samtools-1.8/samtools sort - -o twilhelma_2014_vs_scaffolds_v1.sorted.bam

~/samtools-1.8/samtools view twilhelma_2014_vs_scaffolds_v1.sorted.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > twilhelma_2014_vs_scaffolds_v1.sorted.hits_from_bam.txt

~/git/lavaLampPlot/hits_to_coverage.py -f twilhelma_scaffolds_v1.fasta -b twilhelma_2014_vs_scaffolds_v1.sorted.hits_from_bam.txt > twilhelma_2014_vs_scaffolds_v1.coverage.tab

### v2.0 DNA MAPPING FOR COVERAGE

cd ../twilhelma_v2_dovetail

hisat2-build sponge_tethya_wilhelma_26Apr2018_OM745.fasta sponge_tethya_wilhelma_26Apr2018_OM745.fasta

hisat2 -q -x sponge_tethya_wilhelma_26Apr2018_OM745.fasta -1 ../rawdata/twilhelma_2014_dna_SRR2163223_1.fastq.gz -2 ../rawdata/twilhelma_2014_dna_SRR2163223_2.fastq.gz -S twilhelma_2014_dna_vs_dovetail-v2.sam -p 6 --no-spliced-alignment 2> twilhelma_2014_dna_vs_dovetail-v2.log

#~/samtools-1.8/samtools view -bS twilhelma_2014_dna_vs_dovetail.sam > twilhelma_2014_dna_vs_dovetail.bam
#~/samtools-1.8/samtools sort twilhelma_2014_dna_vs_dovetail.bam -o twilhelma_2014_dna_vs_dovetail.sorted.bam

~/samtools-1.8/samtools view -Su twilhelma_2014_dna_vs_dovetail.sam | ~/samtools-1.8/samtools sort - -o twilhelma_2014_dna_vs_dovetail.sorted.bam

~/samtools-1.8/samtools view twilhelma_2014_dna_vs_dovetail.sorted.bam | sortreadsfrombam.py -i - > twilhelma_2014_dna_vs_dovetail_hits_from_bam.txt


bedgraph2histo.py -b twilhelma_dovetail_bedtools_genomecov.bg > twilhelma_dovetail_bedtools_genomecov.histo.tab

~/samtools-1.8/samtools view twilhelma_2014_dna_vs_dovetail.sorted.bam | sortreadsfrombam.py -i - > twilhelma_2014_dna_vs_dovetail_hits_from_bam.txt

hits_to_coverage.py sponge_tethya_wilhelma_26Apr2018_OM745.fasta twilhelma_2014_dna_vs_dovetail_hits_from_bam.txt > twilhelma_2014_dna_vs_dovetail_coverage.tab

Rscript ~/git/lavaLampPlot/contig_gc_coverage.R twilhelma_2014_dna_vs_dovetail-v2_coverage.tab

### WITHOUT ALLOWING SPLICE ALIGNMENT

hisat2 -q -x sponge_tethya_wilhelma_26Apr2018_OM745.fasta -1 ../rawdata/twilhelma_2014_dna_SRR2163223_1.fastq.gz -2 ../rawdata/twilhelma_2014_dna_SRR2163223_2.fastq.gz -S twilhelma_2014_dna_vs_dovetail-v2.sam -p 6 --no-spliced-alignment 2> twilhelma_2014_dna_vs_dovetail-v2.log

~/samtools-1.8/samtools view -Su twilhelma_2014_dna_vs_dovetail-v2.sam | ~/samtools-1.8/samtools sort - -o twilhelma_2014_dna_vs_dovetail-v2.sorted.bam

~/samtools-1.8/samtools view twilhelma_2014_dna_vs_dovetail-v2.sorted.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > twilhelma_2014_dna_vs_dovetail-v2.sorted.hits_from_bam.txt

### GENERATION OF FILTERED SET v2.1

~/git/lavaLampPlot/hits_to_coverage.py -f sponge_tethya_wilhelma_26Apr2018_OM745.fasta -b twilhelma_2014_dna_vs_dovetail_hits_from_bam.txt -A 20 -W 50 > twilhelma_2014_dna_vs_dovetail_filtered_coverage.tab

Rscript ~/git/lavaLampPlot/contig_gc_coverage.R twilhelma_2014_dna_vs_dovetail_filtered_coverage.tab

cut -f 1 twilhelma_2014_dna_vs_dovetail_filtered_coverage.tab > twilhelma_dovetail_filtered_names

getAinB.py twilhelma_dovetail_filtered_names sponge_tethya_wilhelma_26Apr2018_OM745.fasta > tethya_wilhelma_dovetail.v2.1.fasta

### RNA-SEQ MAPPING FOR TXOME

#v2.0
hisat2 -q -x sponge_tethya_wilhelma_26Apr2018_OM745.fasta -1 ../rawdata/twilhelma_rnaseq_SRR4255675_1.fastq.gz -2 ../rawdata/twilhelma_rnaseq_SRR4255675_2.fastq.gz --rna-strandness RF -S twilhelma_rnaseq_vs_dovetail.sam -p 6 --dta 2> twilhelma_rnaseq_vs_dovetail.log

~/samtools-1.8/samtools view -Su twilhelma_rnaseq_vs_dovetail.sam | ~/samtools-1.8/samtools sort - -o twilhelma_rnaseq_vs_dovetail.sorted.bam

~/stringtie-1.3.4d.Linux_x86_64/stringtie  twilhelma_rnaseq_vs_dovetail.sorted.bam -v -o twilhelma_rnaseq_vs_dovetail.stringtie.gtf -l Tw_v2 -p 2


#v2.1
hisat2 -q -x tethya_wilhelma_dovetail.v2.1.fasta -1 ../rawdata/twilhelma_rnaseq_SRR4255675_1.fastq.gz -2 ../rawdata/twilhelma_rnaseq_SRR4255675_2.fastq.gz --rna-strandness RF -S twilhelma_rnaseq_vs_dovetail.v2.1.sam -p 6 --dta 2> twilhelma_rnaseq_vs_dovetail.v2.1.log

~/samtools-1.8/samtools view -Su twilhelma_rnaseq_vs_dovetail.v2.1.sam | ~/samtools-1.8/samtools sort - -o twilhelma_rnaseq_vs_dovetail.v2.1.sorted.bam

~/stringtie-1.3.4d.Linux_x86_64/stringtie  twilhelma_rnaseq_vs_dovetail.v2.1.sorted.bam -v -o twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.gtf -l Tw_v2 -p 2

~/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.gtf tethya_wilhelma_dovetail.v2.1.fasta > twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta

~/TransDecoder-3.0.1/TransDecoder.LongOrfs -t twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta

~/TransDecoder-3.0.1/TransDecoder.Predict  -t twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta

~/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.gtf > twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.gff

~/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta.transdecoder.gff3 twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.gff twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta > twilhelma_rnaseq_vs_dovetail.v2.1.stringtie.fasta.transdecoder.genome.gff3

```
