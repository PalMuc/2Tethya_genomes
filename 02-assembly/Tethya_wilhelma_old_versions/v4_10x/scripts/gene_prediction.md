# Trimming
Trimming RNA paired end reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```bash 
trimmomatic PE Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq.gz \
Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq.gz \
output_forward_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq.gz \
output_forward_unpaired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq.gz \
output_reverse_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq.gz \
output_reverse_unpaired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq.gz \
HEADCROP:15 
```
## Mapping reads 
RNA reads were mapped to the reference genome using [hisat2](http://daehwankimlab.github.io/hisat2/). Alignments were converted to *bam-files using [samtools](https://github.com/samtools/samtools)
```bash
#!/bin/bash
#
#SBATCH --job-name=hisat2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --time=20:00:00
#SBATCH --mem=55000


hisat2 --avoid-pseudogene \
--novel-splicesite-outfile \
--new-summary \
-x Tethya_wilhelma_V4_below_50_percent_GC \
-1 output_forward_paired_Tethya_RNA-Seq_Fastq1_TAGCTT_lane2.fastq.gz \
-2 output_reverse_paired_Tethya_RNA-Seq_Fastq2_TAGCTT_lane2.fastq.gz \
-S RNA_aligned_T_wilhelma_V4.sam


samtools flagstat  RNA_aligned_T_wilhelma_V4.sam >flagstat.txt

samtools view -F 4 -T Tethya_wilhelma_V4_below_50_percent_GC.fasta -Sb RNA_aligned_T_wilhelma_V4.sam | samtools sort -@ 56 > Twilhelma_V4_mapped_RNA.bam

bam2hints --in=Twilhelma_V4_mapped_RNA.bam --out=hints.gff

```
## Statistics 
```bash
samtools flagstat Twilhelma_V4_mapped_RNA.bam >flagstat.txt
```
```bash
cat flagstat.txt 
216828847 + 0 in total (QC-passed reads + QC-failed reads)
15377273 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
174641042 + 0 mapped (80.54% : N/A)
201451574 + 0 paired in sequencing
100725787 + 0 read1
100725787 + 0 read2
146220276 + 0 properly paired (72.58% : N/A)
149692764 + 0 with itself and mate mapped
9571005 + 0 singletons (4.75% : N/A)
490122 + 0 with mate mapped to a different chr
396896 + 0 with mate mapped to a different chr (mapQ>=5)
```
```perl
#!/usr/bin/perl -w 


open(DAT,"<Tethya_wilhelma_V4_below_50_percent_GC.fasta");
open(OUT,">>modif_Tethya_wilhelma_V4_below_50_percent_GC.fa");
while(<DAT>)
	{
		my$line=$_;
		chomp$line;
		print OUT $line;
		print OUT "\n";
	}

```
## Generating trainingset

Braker2 was used to run Genmark for training-sets.

```bash 
#!/bin/bash
#
#SBATCH --job-name=braker
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --mem=55000


braker.pl --genome modif_Tethya_wilhelma_V4_below_50_percent_GC.fa --bam Twilhelma_V4_mapped_RNA.bam --softmasking --esmode --cores 28


``` 
## BRAKER output
```bash
#**********************************************************************************
#                               BRAKER CONFIGURATION                               
#**********************************************************************************
# BRAKER CALL: /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/braker.pl --genome modif_Tethya_wilhelma_V4_below_50_percent_GC.fa --bam Twilhelma_V4_mapped_RNA.bam --softmasking --cores 28
# Mon May 17 18:38:59 2021: braker.pl version 2.1.5
# Mon May 17 18:38:59 2021: Configuring of BRAKER for using external tools...
# Mon May 17 18:38:59 2021: Found environment variable $AUGUSTUS_CONFIG_PATH. Setting $AUGUSTUS_CONFIG_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/config/
# Mon May 17 18:38:59 2021: Found environment variable $AUGUSTUS_BIN_PATH. Setting $AUGUSTUS_BIN_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/
# Mon May 17 18:38:59 2021: Found environment variable $AUGUSTUS_SCRIPTS_PATH. Setting $AUGUSTUS_SCRIPTS_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/
# Mon May 17 18:38:59 2021: Found environment variable $GENEMARK_PATH. Setting $GENEMARK_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/gmes_linux_64/
# Mon May 17 18:38:59 2021: Did not find environment variable $BAMTOOLS_PATH (either variable does not exist, or the path given in variable does not exist). Will try to set this variable in a different way, later.
# Mon May 17 18:38:59 2021: Trying to guess $BAMTOOLS_BIN_PATH from location of bamtools executable that is available in your $PATH
# Mon May 17 18:38:59 2021: Setting $BAMTOOLS_BIN_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin
# Mon May 17 18:38:59 2021: Did not find environment variable $SAMTOOLS_PATH  (either variable does not exist, or the path given in variable doesnot exist). Will try to set this variable in a different way, later.
# Mon May 17 18:38:59 2021: Trying to guess $SAMTOOLS_PATH from location of samtools executable in your $PATH
# Mon May 17 18:38:59 2021: Setting $SAMTOOLS_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin
# Mon May 17 18:38:59 2021: Trying to guess $DIAMOND_PATH from location of diamond executable that is available in your $PATH
# Mon May 17 18:38:59 2021: Setting $DIAMOND_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin
# Mon May 17 18:38:59 2021: Found environment variable $PYTHON3_PATH. Setting $PYTHON3_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin
# Mon May 17 18:38:59 2021: Did not find environment variable $CDBTOOLS_PATH
# Mon May 17 18:38:59 2021: Trying to guess $CDBTOOLS_PATH from location of cdbfasta executable that is available in your $PATH
# Mon May 17 18:38:59 2021: Setting $CDBTOOLS_PATH to /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin
#*********
# WARNING: Your system does not have 28 cores available, only 1
. Braker will use the 1
  available instead of the chosen 28.
#*********
#*********
# IMPORTANT INFORMATION: no species for identifying the AUGUSTUS  parameter set that will arise from this BRAKER run was set. BRAKER will create an AUGUSTUS parameter set with name Sp_1. This parameter set can be used for future BRAKER/AUGUSTUS prediction runs for the same species. It is usually not necessary to retrain AUGUSTUS with novel extrinsic data if a high quality parameter set already exists.
#*********
#**********************************************************************************
#                               CREATING DIRECTORY STRUCTURE                       
#**********************************************************************************
# Mon May 17 18:39:01 2021: create working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker.
mkdir /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker
# Mon May 17 18:39:01 2021: creating file that contains citations for this BRAKER run at /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/what-to-cite.txt...
# Mon May 17 18:39:01 2021: create working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET.
mkdir /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET
# Mon May 17 18:39:01 2021: create working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/species
mkdir /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/species
# Mon May 17 18:39:01 2021: create working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors
mkdir /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors
# Mon May 17 18:39:02 2021: changing into working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker
cd /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker
# Mon May 17 18:39:02 2021: Creating parameter template files for AUGUSTUS with new_species.pl
# Mon May 17 18:39:02 2021: new_species.pl will create parameter files for species Sp_1 in /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/config//species/Sp_1
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/new_species.pl --species=Sp_1 --AUGUSTUS_CONFIG_PATH=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/config/ 1> /dev/null 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/new_species.stderr
# Mon May 17 18:39:03 2021: check_fasta_headers(): Checking fasta headers of file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/modif_Tethya_wilhelma_V4_below_50_percent_GC.fa
#**********************************************************************************
#                               PROCESSING HINTS                                   
#**********************************************************************************
# Mon May 17 18:39:05 2021: Converting bam files to hints
# Mon May 17 18:39:06 2021: Checking bam headers
# Mon May 17 18:39:06 2021: create header file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/Twilhelma_V4_mapped_RNA_header.sam
/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/bamtools header -in /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/Twilhelma_V4_mapped_RNA.bam > /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/Twilhelma_V4_mapped_RNA_header.sam
# Mon May 17 18:39:06 2021: Deleting SAM header file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/Twilhelma_V4_mapped_RNA_header.sam (will not be needed from here on)
# Mon May 17 18:39:06 2021: Deleting /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/Twilhelma_V4_mapped_RNA_new_header.sam
# Mon May 17 18:39:06 2021: make hints from BAM file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/Twilhelma_V4_mapped_RNA.bam
/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin//bam2hints --intronsonly --in=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/Twilhelma_V4_mapped_RNA.bam --out=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/bam2hints.temp.0.gff 1> /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/bam2hints.0.stdout 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/bam2hints.0.stderr
# Mon May 17 19:25:29 2021: add hints from BAM file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/Twilhelma_V4_mapped_RNA.bam to hints file
cat /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/bam2hints.temp.0.gff >>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.temp.gff
# Mon May 17 19:25:29 2021: Checking for hints of src=C and with grp tags that should not be joined according to multiplicity
# Mon May 17 19:25:29 2021: Joining hints that are identical (& from the same source) into multiplicity hints (input file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/tmp_merge_hints.gff)
# Mon May 17 19:25:29 2021: sort hints of type rnaseq
cat /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/tmp_merge_hints.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 >/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hints.rnaseq.temp.sort.gff
# Mon May 17 19:26:00 2021: join multiple hints
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/join_mult_hints.pl </dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hints.rnaseq.temp.sort.gff >/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/tmp_merge_hints.gff 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/join_mult_hints.rnaseq.stderr
mv /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/tmp_merge_hints.gff /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.temp.gff
# Mon May 17 19:26:14 2021: filter introns, find strand and change score to 'mult' entry
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/filterIntronsFindStrand.pl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genome.fa /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.temp.gff --score 1>>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.gff 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/filterIntronsFindStrand.stderr
# Mon May 17 19:26:27 2021: rm /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.temp.gff
# Mon May 17 19:26:27 2021: Preparing hints for running GeneMark
# Mon May 17 19:26:27 2021: Filtering intron hints for GeneMark from /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.gff...
cat /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff.rnaseq | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/join_mult_hints.pl > /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff.rnaseq.tmp
mv /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff.rnaseq.tmp /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff
#**********************************************************************************
#                              RUNNING GENEMARK-EX                                 
#**********************************************************************************
# Mon May 17 19:27:06 2021: Checking whether file /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff contains sufficient multiplicity information...
# Mon May 17 19:27:09 2021: Executing GeneMark-ET
# Mon May 17 19:27:09 2021: changing into GeneMark-ET directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET
cd /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET
# Mon May 17 19:27:09 2021: Executing gmes_petap.pl
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/program/gmes_linux_64//gmes_petap.pl --verbose --sequence=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genome.fa --ET=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/genemark_hintsfile.gff --et_score 10 --max_intergenic 50000 --cores=28 --gc_donor 0.001 --soft_mask 1000 1>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET.stdout 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/GeneMark-ET.stderr
# Mon May 17 21:47:35 2021: change to working directory /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker
cd /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker
# Mon May 17 21:47:35 2021: Filtering output of GeneMark for generating training genes for AUGUSTUS
# Mon May 17 21:47:35 2021: Checking whether hintsfile contains single exon CDSpart hints or start/stop hints
# Mon May 17 21:47:36 2021: filtering GeneMark genes by intron hints
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/filterGenemark.pl --genemark=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET/genemark.gtf --hints=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/hintsfile.gff 1>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/filterGenemark.stdout 2>/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/filterGenemark.stderr
#Mon May 17 21:48:04 2021: downsampling good genemark genes according to poisson distribution with Lambda 2:
perl /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/miniconda3/envs/BRAKER2/bin/downsample_traingenes.pl --in_gtf=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET/genemark.f.good.gtf --out_gtf=/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/GeneMark-ET/genemark.d.gtf --lambda=2 1> /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/downsample_traingenes.log 2> /dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker/errors/downsample_traingenes.err
```
## Training AUGUSTUS
The training-set was used to train [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) 
```bash
(busco) di52zuc@cm2login2:/dss/dssfs02/lwp-dss-0001/pn69xe/pn69xe-dss-0000/di52zuc/hisat2mapping/braker2> cat train.sh 
#!/bin/bash
#
#SBATCH --job-name=augustus_train
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --mem=55000



new_species.pl --species=Twilhelma_V4

etraining --species=Twilhelma_V4 train.gb

augustus --species=Twilhelma_V4 train.gb

optimize_augustus.pl --species=Twilhelma_V4 train.gb

etraining --species=Twilhelma_V4 train.gb 
```

## Protein prediction 

```bash
#!/bin/bash
#
#SBATCH --job-name=augustus_prot_pred
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp3
#SBATCH --time=20:00:00
#SBATCH --mem=55000

augustus --strand=both \
--codingseq=on \
--genemodel=partial \
--species=Twilhelma_V4 \
--outfile=predicted_genes_coding_sequences.gff \
genome.fa
```
## Genes per scaffold
| scaffold ID 	| genes 	|
|-	|-	|
| P_RNA_scaffold_29 	| 6079 	|
| P_RNA_scaffold_27 	| 4340 	|
| P_RNA_scaffold_10 	| 2950 	|
| scaffold6 	| 1449 	|
| P_RNA_scaffold_2 	| 1427 	|
| P_RNA_scaffold_44 	| 1398 	|
| P_RNA_scaffold_36 	| 1301 	|
| P_RNA_scaffold_24 	| 1287 	|
| P_RNA_scaffold_43 	| 970 	|
| P_RNA_scaffold_5 	| 655 	|
| P_RNA_scaffold_38 	| 567 	|
| P_RNA_scaffold_13 	| 379 	|
| P_RNA_scaffold_8 	| 300 	|
| scaffold20 	| 185 	|
| P_RNA_scaffold_6 	| 173 	|
| P_RNA_scaffold_28 	| 161 	|
| scaffold24 	| 141 	|
| P_RNA_scaffold_18 	| 134 	|
| P_RNA_scaffold_33 	| 96 	|
| P_RNA_scaffold_7 	| 82 	|
| P_RNA_scaffold_21 	| 73 	|
| scaffold59 	| 72 	|
| scaffold31 	| 66 	|
| P_RNA_scaffold_11 	| 58 	|
| scaffold34 	| 57 	|
| scaffold41 	| 56 	|
| scaffold43 	| 54 	|
| scaffold46 	| 49 	|
| P_RNA_scaffold_19 	| 49 	|
| scaffold45 	| 48 	|
| scaffold39 	| 47 	|
| P_RNA_scaffold_48 	| 46 	|
| scaffold48 	| 45 	|
| P_RNA_scaffold_1 	| 44 	|
| scaffold54 	| 44 	|
| scaffold52 	| 44 	|
| P_RNA_scaffold_16 	| 44 	|
| P_RNA_scaffold_22 	| 42 	|
| P_RNA_scaffold_39 	| 42 	|
| scaffold47 	| 39 	|
| P_RNA_scaffold_35 	| 37 	|
| scaffold92 	| 35 	|
| P_RNA_scaffold_30 	| 34 	|
| P_RNA_scaffold_12 	| 32 	|
| scaffold68 	| 32 	|
| scaffold65 	| 32 	|
| scaffold80 	| 31 	|
| scaffold127 	| 31 	|
| P_RNA_scaffold_26 	| 30 	|
| scaffold60 	| 29 	|
| scaffold57 	| 28 	|
| scaffold86 	| 27 	|
| scaffold82 	| 26 	|
| scaffold64 	| 26 	|
| scaffold97 	| 25 	|
| scaffold110 	| 25 	|
| scaffold53 	| 25 	|
| scaffold107 	| 24 	|
| scaffold137 	| 24 	|
| P_RNA_scaffold_34 	| 24 	|
| scaffold96 	| 23 	|
| P_RNA_scaffold_3 	| 23 	|
| scaffold120 	| 23 	|
| scaffold116 	| 22 	|
| scaffold111 	| 22 	|
| scaffold67 	| 22 	|
| P_RNA_scaffold_41 	| 21 	|
| P_RNA_scaffold_4 	| 20 	|
| scaffold74 	| 20 	|
| P_RNA_scaffold_9 	| 20 	|
| scaffold84 	| 19 	|
| P_RNA_scaffold_46 	| 19 	|
| scaffold89 	| 19 	|
| scaffold134 	| 19 	|
| scaffold99 	| 19 	|
| scaffold85 	| 18 	|
| scaffold154 	| 18 	|
| scaffold69 	| 18 	|
| scaffold208 	| 17 	|
| scaffold147 	| 17 	|
| scaffold156 	| 17 	|
| scaffold114 	| 17 	|
| P_RNA_scaffold_45 	| 17 	|
| scaffold148 	| 16 	|
| scaffold175 	| 16 	|
| scaffold159 	| 16 	|
| scaffold140 	| 16 	|
| scaffold106 	| 16 	|
| P_RNA_scaffold_31 	| 16 	|
| scaffold91 	| 16 	|
| scaffold155 	| 15 	|
| scaffold153 	| 15 	|
| scaffold115 	| 15 	|
| scaffold126 	| 15 	|
| scaffold165 	| 14 	|
| scaffold171 	| 14 	|
| P_RNA_scaffold_17 	| 14 	|
| scaffold170 	| 14 	|
| scaffold157 	| 14 	|
| scaffold160 	| 14 	|
| P_RNA_scaffold_37 	| 14 	|
| scaffold142 	| 14 	|
| scaffold180 	| 13 	|
| scaffold179 	| 13 	|
| scaffold258 	| 12 	|
| scaffold145 	| 12 	|
| scaffold150 	| 12 	|
| scaffold131 	| 12 	|
| scaffold194 	| 11 	|
| scaffold213 	| 10 	|
| scaffold206 	| 10 	|
| scaffold215 	| 10 	|
| scaffold217 	| 10 	|
| scaffold209 	| 10 	|
| scaffold244 	| 10 	|
| scaffold317 	| 10 	|
| scaffold141 	| 10 	|
| scaffold224 	| 10 	|
| scaffold149 	| 10 	|
| scaffold204 	| 9 	|
| scaffold218 	| 9 	|
| scaffold190 	| 9 	|
| scaffold214 	| 9 	|
| scaffold184 	| 9 	|
| scaffold221 	| 9 	|
| P_RNA_scaffold_42 	| 9 	|
| scaffold187 	| 9 	|
| scaffold201 	| 9 	|
| scaffold167 	| 9 	|
| scaffold291 	| 9 	|
| scaffold399 	| 8 	|
| scaffold188 	| 8 	|
| scaffold301 	| 8 	|
| scaffold228 	| 8 	|
| P_RNA_scaffold_20 	| 8 	|
| scaffold196 	| 8 	|
| scaffold259 	| 8 	|
| scaffold138 	| 8 	|
| scaffold212 	| 8 	|
| scaffold243 	| 8 	|
| scaffold370 	| 7 	|
| scaffold289 	| 7 	|
| scaffold225 	| 7 	|
| scaffold268 	| 7 	|
| scaffold283 	| 7 	|
| scaffold251 	| 7 	|
| scaffold161 	| 7 	|
| scaffold223 	| 7 	|
| scaffold197 	| 7 	|
| scaffold146 	| 7 	|
| scaffold192 	| 7 	|
| scaffold280 	| 6 	|
| scaffold282 	| 6 	|
| scaffold271 	| 6 	|
| scaffold284 	| 6 	|
| scaffold256 	| 6 	|
| scaffold310 	| 6 	|
| scaffold285 	| 6 	|
| scaffold152 	| 6 	|
| scaffold347 	| 5 	|
| scaffold290 	| 5 	|
| scaffold315 	| 5 	|
| scaffold337 	| 5 	|
| scaffold446 	| 5 	|
| scaffold331 	| 5 	|
| scaffold295 	| 5 	|
| scaffold210 	| 5 	|
| scaffold236 	| 5 	|
| scaffold287 	| 5 	|
| scaffold237 	| 5 	|
| scaffold231 	| 5 	|
| scaffold202 	| 5 	|
| scaffold299 	| 5 	|
| scaffold300 	| 5 	|
| scaffold383 	| 4 	|
| scaffold320 	| 4 	|
| scaffold385 	| 4 	|
| scaffold469 	| 4 	|
| scaffold278 	| 4 	|
| scaffold330 	| 4 	|
| scaffold267 	| 4 	|
| P_RNA_scaffold_40 	| 4 	|
| scaffold269 	| 4 	|
| scaffold325 	| 4 	|
| scaffold238 	| 4 	|
| P_RNA_scaffold_47 	| 4 	|
| scaffold270 	| 4 	|
| scaffold306 	| 3 	|
| scaffold507 	| 3 	|
| scaffold438 	| 3 	|
| scaffold472 	| 3 	|
| scaffold452 	| 3 	|
| scaffold424 	| 3 	|
| scaffold454 	| 3 	|
| scaffold342 	| 3 	|
| scaffold420 	| 3 	|
| scaffold323 	| 3 	|
| scaffold274 	| 3 	|
| scaffold384 	| 3 	|
| scaffold527 	| 3 	|
| scaffold433 	| 3 	|
| scaffold494 	| 3 	|
| scaffold448 	| 3 	|
| scaffold462 	| 3 	|
| scaffold334 	| 3 	|
| scaffold288 	| 3 	|
| scaffold305 	| 3 	|
| scaffold344 	| 3 	|
| scaffold613 	| 3 	|
| scaffold923 	| 2 	|
| scaffold628 	| 2 	|
| scaffold616 	| 2 	|
| scaffold487 	| 2 	|
| scaffold561 	| 2 	|
| scaffold352 	| 2 	|
| scaffold828 	| 2 	|
| scaffold414 	| 2 	|
| scaffold706 	| 2 	|
| scaffold631 	| 2 	|
| scaffold373 	| 2 	|
| scaffold568 	| 2 	|
| scaffold549 	| 2 	|
| scaffold801 	| 2 	|
| scaffold573 	| 2 	|
| scaffold513 	| 2 	|
| scaffold551 	| 2 	|
| scaffold356 	| 2 	|
| scaffold451 	| 2 	|
| scaffold614 	| 2 	|
| scaffold576 	| 2 	|
| scaffold508 	| 2 	|
| scaffold540 	| 2 	|
| scaffold546 	| 2 	|
| scaffold566 	| 2 	|
| scaffold419 	| 2 	|
| scaffold697 	| 2 	|
| scaffold434 	| 2 	|
| scaffold562 	| 2 	|
| scaffold250 	| 2 	|
| scaffold313 	| 2 	|
| scaffold735 	| 2 	|
| scaffold538 	| 2 	|
| P_RNA_scaffold_25 	| 2 	|
| scaffold405 	| 2 	|
| scaffold473 	| 2 	|
| scaffold650 	| 2 	|
| scaffold575 	| 2 	|
| scaffold544 	| 2 	|
| scaffold413 	| 2 	|
| scaffold316 	| 2 	|
| scaffold346 	| 2 	|
| scaffold542 	| 2 	|
| P_RNA_scaffold_23 	| 2 	|
| scaffold222 	| 2 	|
| scaffold491 	| 2 	|
| scaffold595 	| 2 	|
| scaffold474 	| 2 	|
| scaffold496 	| 2 	|
| scaffold353 	| 2 	|
| scaffold457 	| 2 	|
| scaffold391 	| 2 	|
| scaffold681 	| 2 	|
| scaffold633 	| 2 	|
| scaffold583 	| 2 	|
| scaffold671 	| 2 	|
| scaffold460 	| 2 	|
| scaffold791 	| 2 	|
| scaffold874 	| 1 	|
| scaffold928 	| 1 	|
| scaffold798 	| 1 	|
| scaffold483 	| 1 	|
| scaffold788 	| 1 	|
| scaffold630 	| 1 	|
| scaffold759 	| 1 	|
| scaffold690 	| 1 	|
| scaffold800 	| 1 	|
| scaffold858 	| 1 	|
| scaffold803 	| 1 	|
| scaffold478 	| 1 	|
| scaffold955 	| 1 	|
| scaffold668 	| 1 	|
| scaffold557 	| 1 	|
| scaffold696 	| 1 	|
| scaffold710 	| 1 	|
| scaffold622 	| 1 	|
| scaffold773 	| 1 	|
| scaffold634 	| 1 	|
| scaffold941 	| 1 	|
| scaffold809 	| 1 	|
| scaffold902 	| 1 	|
| scaffold642 	| 1 	|
| scaffold812 	| 1 	|
| scaffold962 	| 1 	|
| scaffold813 	| 1 	|
| scaffold652 	| 1 	|
| scaffold679 	| 1 	|
| scaffold484 	| 1 	|
| scaffold826 	| 1 	|
| scaffold914 	| 1 	|
| scaffold911 	| 1 	|
| scaffold934 	| 1 	|
| scaffold605 	| 1 	|
| scaffold898 	| 1 	|
| scaffold865 	| 1 	|
| scaffold965 	| 1 	|
| scaffold721 	| 1 	|
| scaffold688 	| 1 	|
| scaffold714 	| 1 	|
| scaffold931 	| 1 	|
| scaffold856 	| 1 	|
| scaffold689 	| 1 	|
| scaffold669 	| 1 	|
| scaffold730 	| 1 	|
| scaffold783 	| 1 	|
| scaffold465 	| 1 	|
| scaffold511 	| 1 	|
| scaffold954 	| 1 	|
| scaffold427 	| 1 	|
| scaffold637 	| 1 	|
| scaffold732 	| 1 	|
| scaffold657 	| 1 	|
| scaffold933 	| 1 	|
| scaffold761 	| 1 	|
| scaffold476 	| 1 	|
| scaffold619 	| 1 	|
| scaffold559 	| 1 	|
| scaffold638 	| 1 	|
| scaffold820 	| 1 	|
| scaffold685 	| 1 	|
| scaffold428 	| 1 	|
| scaffold754 	| 1 	|
| scaffold646 	| 1 	|
| scaffold845 	| 1 	|
| scaffold700 	| 1 	|
| scaffold321 	| 1 	|
| scaffold765 	| 1 	|
| scaffold779 	| 1 	|
| scaffold794 	| 1 	|
| scaffold797 	| 1 	|
| scaffold379 	| 1 	|
| scaffold517 	| 1 	|
| scaffold623 	| 1 	|
| scaffold917 	| 1 	|
| scaffold403 	| 1 	|
| scaffold663 	| 1 	|
| scaffold456 	| 1 	|
| scaffold624 	| 1 	|
| scaffold882 	| 1 	|
| scaffold664 	| 1 	|
| scaffold577 	| 1 	|
| scaffold695 	| 1 	|
| scaffold459 	| 1 	|
| scaffold587 	| 1 	|
| scaffold831 	| 1 	|
| scaffold885 	| 1 	|
| scaffold640 	| 1 	|
| scaffold503 	| 1 	|
| scaffold804 	| 1 	|
| scaffold441 	| 1 	|
| scaffold589 	| 1 	|
| scaffold742 	| 1 	|
| scaffold698 	| 1 	|
| scaffold945 	| 1 	|
| scaffold751 	| 1 	|
| scaffold564 	| 1 	|
| scaffold939 	| 1 	|
| scaffold903 	| 1 	|
| scaffold753 	| 1 	|
| scaffold659 	| 1 	|
| scaffold308 	| 1 	|
| scaffold745 	| 1 	|
| scaffold960 	| 1 	|
| scaffold770 	| 1 	|
| scaffold862 	| 1 	|
| scaffold716 	| 1 	|
| scaffold760 	| 1 	|
| scaffold835 	| 1 	|
| scaffold541 	| 1 	|
| scaffold429 	| 1 	|
| scaffold755 	| 1 	|
| scaffold550 	| 1 	|
| scaffold864 	| 1 	|
| scaffold611 	| 1 	|
| scaffold924 	| 1 	|
| scaffold588 	| 1 	|
| scaffold603 	| 1 	|
| scaffold529 	| 1 	|
| scaffold617 	| 1 	|
| scaffold506 	| 1 	|
| scaffold711 	| 1 	|
| scaffold302 	| 1 	|
| scaffold781 	| 1 	|
| scaffold596 	| 1 	|
| scaffold660 	| 1 	|
| scaffold739 	| 1 	|
| scaffold883 	| 1 	|
| scaffold707 	| 1 	|
| scaffold784 	| 1 	|
| scaffold769 	| 1 	|
| scaffold606 	| 1 	|
| scaffold612 	| 1 	|
| scaffold591 	| 1 	|
| scaffold398 	| 1 	|
| scaffold913 	| 1 	|
| scaffold785 	| 1 	|
| scaffold329 	| 1 	|
| scaffold719 	| 1 	|
| scaffold504 	| 1 	|
| scaffold548 	| 1 	|
| scaffold808 	| 1 	|
| scaffold921 	| 1 	|
| scaffold774 	| 1 	|
| scaffold524 	| 1 	|
| scaffold627 	| 1 	|
| scaffold678 	| 1 	|
| scaffold814 	| 1 	|
| scaffold836 	| 1 	|
| scaffold802 	| 1 	|
| scaffold470 	| 1 	|
| scaffold901 	| 1 	|
| scaffold786 	| 1 	|
| scaffold869 	| 1 	|
| scaffold958 	| 1 	|
| scaffold666 	| 1 	|
| scaffold449 	| 1 	|
| scaffold667 	| 1 	|
| P_RNA_scaffold_32 	| 1 	|
| scaffold363 	| 1 	|
| scaffold672 	| 1 	|
| scaffold793 	| 1 	|
| scaffold846 	| 1 	|
| scaffold380 	| 1 	|
| scaffold823 	| 0 	|
| scaffold893 	| 0 	|
| scaffold942 	| 0 	|
| scaffold909 	| 0 	|
| scaffold891 	| 0 	|
| scaffold810 	| 0 	|
| scaffold870 	| 0 	|
| scaffold872 	| 0 	|
| scaffold957 	| 0 	|
| scaffold946 	| 0 	|
| scaffold861 	| 0 	|
| scaffold847 	| 0 	|
| scaffold848 	| 0 	|
| scaffold834 	| 0 	|
| scaffold539 	| 0 	|
| scaffold927 	| 0 	|
| scaffold964 	| 0 	|
| scaffold839 	| 0 	|
| scaffold677 	| 0 	|
| scaffold776 	| 0 	|
| scaffold737 	| 0 	|
| scaffold961 	| 0 	|
| scaffold886 	| 0 	|
| scaffold837 	| 0 	|
| scaffold853 	| 0 	|
| scaffold694 	| 0 	|
| scaffold593 	| 0 	|
| scaffold670 	| 0 	|
| scaffold560 	| 0 	|
| scaffold653 	| 0 	|
| scaffold674 	| 0 	|
| scaffold686 	| 0 	|
| scaffold553 	| 0 	|
| scaffold523 	| 0 	|
| scaffold599 	| 0 	|
| scaffold827 	| 0 	|
| scaffold704 	| 0 	|
| scaffold376 	| 0 	|
| scaffold767 	| 0 	|
| scaffold873 	| 0 	|
| scaffold740 	| 0 	|
| scaffold881 	| 0 	|
| scaffold775 	| 0 	|
| scaffold925 	| 0 	|
| scaffold904 	| 0 	|
| scaffold790 	| 0 	|
| scaffold818 	| 0 	|
| scaffold665 	| 0 	|
| scaffold948 	| 0 	|
| scaffold485 	| 0 	|
| scaffold377 	| 0 	|
| scaffold947 	| 0 	|
| scaffold608 	| 0 	|
| scaffold866 	| 0 	|
| scaffold632 	| 0 	|
| scaffold772 	| 0 	|
| scaffold876 	| 0 	|
| scaffold887 	| 0 	|
| scaffold838 	| 0 	|
| scaffold708 	| 0 	|
| scaffold768 	| 0 	|
| scaffold789 	| 0 	|
| scaffold703 	| 0 	|
| scaffold944 	| 0 	|
| scaffold877 	| 0 	|
| scaffold930 	| 0 	|
| scaffold675 	| 0 	|
| scaffold725 	| 0 	|
| scaffold943 	| 0 	|
| scaffold918 	| 0 	|
| scaffold932 	| 0 	|
| scaffold897 	| 0 	|
| scaffold729 	| 0 	|
| scaffold763 	| 0 	|
| scaffold416 	| 0 	|
| scaffold879 	| 0 	|
| scaffold822 	| 0 	|
| scaffold578 	| 0 	|
| scaffold649 	| 0 	|
| scaffold691 	| 0 	|
| scaffold796 	| 0 	|
| scaffold750 	| 0 	|
| scaffold626 	| 0 	|
| scaffold722 	| 0 	|
| scaffold458 	| 0 	|
| scaffold824 	| 0 	|
| scaffold787 	| 0 	|
| scaffold830 	| 0 	|
| scaffold821 	| 0 	|
| scaffold895 	| 0 	|
| scaffold684 	| 0 	|
| scaffold743 	| 0 	|
| scaffold635 	| 0 	|
| scaffold764 	| 0 	|
| scaffold833 	| 0 	|
| scaffold890 	| 0 	|
| scaffold855 	| 0 	|
| scaffold586 	| 0 	|
| scaffold840 	| 0 	|
| scaffold579 	| 0 	|
| scaffold816 	| 0 	|
| scaffold658 	| 0 	|
| scaffold899 	| 0 	|
| scaffold766 	| 0 	|
| scaffold749 	| 0 	|
| scaffold481 	| 0 	|
| scaffold959 	| 0 	|
| scaffold852 	| 0 	|
| scaffold938 	| 0 	|
| scaffold699 	| 0 	|
| scaffold953 	| 0 	|
| scaffold594 	| 0 	|
| scaffold907 	| 0 	|
| scaffold734 	| 0 	|
| scaffold849 	| 0 	|
| scaffold718 	| 0 	|
| scaffold829 	| 0 	|
| scaffold951 	| 0 	|
| scaffold705 	| 0 	|
| scaffold807 	| 0 	|
| scaffold929 	| 0 	|
| scaffold811 	| 0 	|
| scaffold956 	| 0 	|
| scaffold963 	| 0 	|
| scaffold900 	| 0 	|
| scaffold604 	| 0 	|
| scaffold641 	| 0 	|
| scaffold949 	| 0 	|
| scaffold281 	| 0 	|
| scaffold748 	| 0 	|
| scaffold582 	| 0 	|
| scaffold967 	| 0 	|
| scaffold851 	| 0 	|
| scaffold795 	| 0 	|
| scaffold842 	| 0 	|
| scaffold910 	| 0 	|
| scaffold662 	| 0 	|
| scaffold940 	| 0 	|
| scaffold567 	| 0 	|
| scaffold936 	| 0 	|
| scaffold727 	| 0 	|
| scaffold350 	| 0 	|
| scaffold521 	| 0 	|
| scaffold528 	| 0 	|
| scaffold570 	| 0 	|

