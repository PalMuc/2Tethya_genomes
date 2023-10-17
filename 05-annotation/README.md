# *Tethya wilhelma* annotation steps #

```
~/git/genomeGTFtools/rename_gtf_contigs.py --omit-comments -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g predicted_genes_Tethya_wilhelma_V4_below_50_perc.gff.gz > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff

extract_augustus_features.py -n TwiV4 -c TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.nucl.fasta -p TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.prot.fasta TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff

~/git/genomeGTFtools/misc/augustus_to_gff3.py TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff

~/gffread-0.12.7.Linux_x86_64/gffread -g Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta -w TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff
prottrans.py -n -r TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta | remove_identical_seqs.py - > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta

~/diamond-v2.0.13/diamond-v2.0.13 makedb --in TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta

cat TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.w_decoys.fasta
~/salmon-latest_linux_x86_64/bin/salmon index -t TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.w_decoys.fasta -i TwiV4_AUG_tx_index -d decoys.txt -k 31 -p 4
```

###  *T. wilhelma* v3 annotation ###

```
perl ../program/BRAKER_v2.0/braker.pl --GENEMARK_PATH=/user/gm_et_linux_64/gmes_petap/ \
--AUGUSTUS_CONFIG_PATH=/user/augustus-3.3.1/config/ \
--AUGUSTUS_BIN_PATH=/user/augustus-3.3.1/bin/ \
--AUGUSTUS_SCRIPTS_PATH=/user/augustus-3.3.1/scripts/ \
--BAMTOOLS_PATH=/user/bamtools-master/build/src/toolkit/ \
--useexisting \
--species=Tetya_wilh \
--genome=gaps_T-wilhelma_after_HM_and_SSPACE.scaffolds.fasta --bam=T_wilhelma_mapped_RNA.bam 

augustus-3.3.1/bin/augustus --strand=both --codingseq=on --genemodel=partial --AUGUSTUS_CONFIG_PATH=/user/augustus-3.3.1/config --species=Tethya_wilh --outfile=predicted_genes_gaps_closed_scaffolded_10_x_raw_coding_sequences.gff gaps_closed_scaffolded_10_x_raw.fasta

```

# *Tethya minuta* annotation #
Annotation of the filtered version was made at [WebAUGUSTUS](https://bioinf.uni-greifswald.de/webaugustus/index), with [resuts here](https://bioinf.uni-greifswald.de/webaugustus/training/show/0db58ea18774e3e101877f3fe1bf0011).

### *T. minuta* first annotation unfiltered #
Trying to make dotplot with [dgenies](https://dgenies.toulouse.inrae.fr/install):

`~/gffread-0.12.7.Linux_x86_64/gffread -g TmiV4_sponge_bin.17.fa -w Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.nucl.fasta Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.gff`

`prottrans.py -n -r Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.nucl.fasta | remove_identical_seqs.py - > Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.prot.fasta`

Mapping was insufficient, trying with proteins:


`~/gffread-0.12.7.Linux_x86_64/gffread -g Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta -w TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff`

`prottrans.py -n -r TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta | remove_identical_seqs.py - > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta`

```
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.prot.fasta -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta -o Tmi_AUG_V4b17_vs_TwiV4_AUG.blastp.tab
#Reported 63580 pairwise alignments, 63580 HSPs.
#7733 queries aligned.
```

### *T. minuta* AUGUSTUS annotation without bacteria ###

Trying scaffold synteny with `V4b17 AUGUSTUS`:

```
~/git/genomeGTFtools/scaffold_synteny.py -b Tmi_AUG_V4b17_vs_TwiV4_AUG.blastp.tab -q Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.gff -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff -f TmiV4_sponge_bin.17.fa -F Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta -l 90 -L 120 --local-positions > Tmi_AUG_V4b17_vs_TwiV4_AUG.blastp.synteny.l.tab
# Parsing genomic contigs TmiV4_sponge_bin.17.fa   17:26:20 2023
# Found 255 contigs   17:26:21 2023
# Sorting contigs by length, keeping up to 90Mbp
# Kept 227 contigs, for 90011345 bases, last contig was 17123bp long   17:26:21 2023
# Parsing genomic contigs Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta   17:26:21 2023
# Found 557 contigs   17:26:21 2023
# Sorting contigs by length, keeping up to 120Mbp
# Kept 76 contigs, for 120045593 bases, last contig was 88337bp long   17:26:21 2023
# Parsing loci from Tethya_minuta_AUGUSTUS_V4.no_comment.sponge_only.gff   17:26:21 2023
# Found 17820 genes   17:26:21 2023
# GFF names parsed as g25173.t1 from g25173.t1
# Parsing loci from TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff   17:26:21 2023
# Found 112344 genes   17:26:22 2023
# GFF names parsed as g28113.t1 from g28113.t1
# Parsing tabular blast output Tmi_AUG_V4b17_vs_TwiV4_AUG.blastp.tab   17:26:22 2023
# Found blast hits for 7710 query sequences, removed 866 hits by evalue   17:26:22 2023
# Removed 0 queries and 0 subjects with 250 or more hits
# Blast names parsed as g25173.t1 from g25173.t1, and g17562.t1 from g17562.t1
# Kept 7710 blast hits
# Determining match positions   17:26:22 2023
# Wrote match positions for 7505 genes
# Wrote target positions for 5170 genes
```

Some scaffolds have no matches due to no AUGUSTUS genes, but have Trinity mapping, including `Scaffolds_1040`

```
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q augustus.hints_utr.aa.fasta -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta -o Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.tab
#Reported 154997 pairwise alignments, 154997 HSPs.
#20556 queries aligned.
```

```
~/git/genomeGTFtools/scaffold_synteny.py -b Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.tab -q augustus.hints_utr.no_comment.gff -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff -f TmiV4_sponge_bin.17.nr.fa -F Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta -E duplicate_scaffolds_to_exclude.names -l 90 -L 120 > Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.synteny.tab
# Reading exclusion list duplicate_scaffolds_to_exclude.names      18:09:55 2023
# Found 11 contigs to exclude      18:09:55 2023
# Parsing genomic contigs TmiV4_sponge_bin.17.nr.fa      18:09:55 2023
# Found 244 contigs      18:09:55 2023
# Sorting contigs by length, keeping up to 90Mbp
# Kept 244 contigs, for 86067122 bases, last contig was 5213bp long      18:09:55 2023
# Parsing genomic contigs Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta      18:09:55 2023
# Found 557 contigs      18:09:56 2023
# Sorting contigs by length, keeping up to 120Mbp
# Kept 76 contigs, for 120045593 bases, last contig was 88337bp long      18:09:56 2023
# Parsing loci from augustus.hints_utr.no_comment.gff      18:09:56 2023
# Found 43394 genes      18:09:57 2023
# GFF names parsed as g22779.t1 from g22779.t1
# 21697 RNA to protein IDs parsed as g22779.t1 from None
# Parsing loci from TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff      18:09:57 2023
# Found 112344 genes      18:09:58 2023
# GFF names parsed as g28113.t1 from g28113.t1
# 1 RNA to protein IDs parsed as g28113.t1 from None
# Parsing tabular blast output Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.tab      18:09:58 2023
# Found blast hits for 20500 query sequences, removed 2262 hits by evalue      18:09:58 2023
# Removed 0 queries and 0 subjects with 250 or more hits
# Blast names parsed as g22779.t1 from g22779.t1, and g23028.t1 from g23028.t1
# Kept 20500 blast hits
# Determining match positions      18:09:58 2023
# Wrote match positions for 19077 genes
# Wrote target positions for 12492 genes

Rscript ~/git/genomeGTFtools/synteny_2d_plot.R Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.synteny.tab Tethya-minuta Tethya-wilhelma 120
```

```
~/git/genomeGTFtools/microsynteny.py -b Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.tab -q augustus.hints_utr.no_comment.gff -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff -E duplicate_scaffolds_to_exclude.names -T > Tmi_V4b17_hintsutr_vs_TwiV4_AUG.microsynteny.T.tab
Rscript ~/git/genomeGTFtools/synteny_block_length_plot.R Tmi_V4b17_hintsutr_vs_TwiV4_AUG.microsynteny.T.tab 0.45
```


# mapping between both genomes #

```
sizecutter.py -Q Tethya_minuta_Trinity.fa
# Counted 151079 sequences
# Total input letters is: 102364741
# Average length is: 677.56
# Median length is: 385
# Shortest sequence is: 183
# Longest sequence is: 15593
# n50 length is: 51183226 and n50 is: 1049.00
# n10 length is: 92128924 and n10 is: 3337.00
```

```
~/minimap2-2.23_x64-linux/minimap2 -a -x splice --secondary=no TmiV4_sponge_bin.17.fa Tethya_minuta_Trinity.fa | ~/samtools-1.14/samtools sort - -o trinity.vs_v4_bin17.bam
[M::mm_idx_gen::3.908*0.83] collected minimizers
[M::mm_idx_gen::4.628*1.17] sorted minimizers
[M::main::4.628*1.17] loaded/built the index for 255 target sequence(s)
[M::mm_mapopt_update::4.950*1.15] mid_occ = 75
[M::mm_idx_stat] kmer size: 15; skip: 5; is_hpc: 0; #seq: 255
[M::mm_idx_stat::5.147*1.15] distinct minimizers: 20269751 (74.64% are singletons); average occurrences: 1.514; average spacing: 2.943; total length: 90289220
[M::worker_pipeline::60.244*2.82] mapped 151079 sequences
[M::main] Version: 2.23-r1111
[M::main] CMD: /home/dummy/minimap2-2.23_x64-linux/minimap2 -a -x splice --secondary=no TmiV4_sponge_bin.17.fa Tethya_minuta_Trinity.fa
[M::main] Real time: 60.274 sec; CPU: 170.017 sec; Peak RSS: 2.456 GB
```

Mapping rate is very low.

```
~/samtools-1.14/samtools flagstat trinity.vs_v4_bin17.bam
152427 + 0 in total (QC-passed reads + QC-failed reads)
151079 + 0 primary
0 + 0 secondary
1348 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
62579 + 0 mapped (41.06% : N/A)
61231 + 0 primary mapped (40.53% : N/A)
```

At least, filtering was likely successful, as very few transcripts map outside of bin 17:

```
~/samtools-1.14/samtools flagstat trinity.vs_v4_all.bam
152439 + 0 in total (QC-passed reads + QC-failed reads)
151079 + 0 primary
0 + 0 secondary
1360 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
62967 + 0 mapped (41.31% : N/A)
61607 + 0 primary mapped (40.78% : N/A)
```

Check assembly against TwiV4

```
~/minimap2-2.23_x64-linux/minimap2 -cx asm5 ../TwiV4.main.fasta TmiV4_sponge_bin.17.fa > TmiV4.17_vs_TwiV4.minimap_asm5.paf
[M::mm_idx_gen::2.299*1.27] collected minimizers
[M::mm_idx_gen::2.595*1.46] sorted minimizers
[M::main::2.595*1.46] loaded/built the index for 557 target sequence(s)
[M::mm_mapopt_update::2.729*1.44] mid_occ = 69
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 557
[M::mm_idx_stat::2.819*1.43] distinct minimizers: 9916314 (89.65% are singletons); average occurrences: 1.283; average spacing: 9.917; total length: 126187550
[M::worker_pipeline::184.324*2.81] mapped 255 sequences
[M::main] Version: 2.23-r1111
[M::main] CMD: /home/wrf/minimap2-2.23_x64-linux/minimap2 -cx asm5 ../TwiV4.main.fasta TmiV4_sponge_bin.17.fa
[M::main] Real time: 184.343 sec; CPU: 518.022 sec; Peak RSS: 2.701 GB
```

Does not work with `asm5` mode, checking with `asm20`:

```
 ~/minimap2-2.23_x64-linux/minimap2 -cx asm20 ../TwiV4.main.fasta TmiV4_sponge_bin.17.fa > TmiV4.17_vs_TwiV4.minimap_asm20.paf
[M::mm_idx_gen::2.832*1.33] collected minimizers
[M::mm_idx_gen::3.349*1.59] sorted minimizers
[M::main::3.349*1.59] loaded/built the index for 557 target sequence(s)
[M::mm_mapopt_update::3.614*1.54] mid_occ = 71
[M::mm_idx_stat] kmer size: 19; skip: 10; is_hpc: 0; #seq: 557
[M::mm_idx_stat::3.799*1.52] distinct minimizers: 17929940 (89.45% are singletons); average occurrences: 1.277; average spacing: 5.510; total length: 126187550
[M::worker_pipeline::35.157*2.53] mapped 255 sequences
[M::main] Version: 2.23-r1111
[M::main] CMD: /home/wrf/minimap2-2.23_x64-linux/minimap2 -cx asm20 ../TwiV4.main.fasta TmiV4_sponge_bin.17.fa
[M::main] Real time: 35.193 sec; CPU: 89.138 sec; Peak RSS: 3.323 GB
```

## generating orthogroups ##

Making the de-duplicated version of the *T.minuta* AUGUSTUS annotation:

```
grep ">" TmiV4_sponge_bin.17.nr.fa | cut -c 2- > TmiV4_sponge_bin.17.nr.names
getgtfscaffolds.py -s TmiV4_sponge_bin.17.nr.names -g augustus.hints_utr.no_comment.gff > augustus.hints_utr.no_comment.nr.gff
# Reading list from TmiV4_sponge_bin.17.nr.names
# Read 244 lines for 244 entries
# Reading features from augustus.hints_utr.no_comment.gff
# Counted 367573 lines, 0 comments
# Wrote 348246 lines
grep mRNA augustus.hints_utr.no_comment.nr.gff | cut -f 9 | cut -d ";" -f 1 | cut -d "=" -f 2 > augustus.hints_utr.aa.nr.names
getAinB.py augustus.hints_utr.aa.nr.names augustus.hints_utr.aa.fasta > augustus.hints_utr.aa.nr.fasta
```

Generate ortholog groups:

```
fastarenamer.py -p "Tmin" -j "|" augustus.hints_utr.aa.nr.fasta > augustus.hints_utr.aa.nr.n.fasta
fastarenamer.py -p "Twil" -j "|" TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.fasta > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.n.fasta
~/diamond-v2.0.13/diamond-v2.0.13 makedb --in TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.n.fasta -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.n.fasta
~/diamond-v2.0.13/diamond-v2.0.13 blastp -q augustus.hints_utr.aa.nr.n.fasta -d TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.n.fasta -o Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.n.tab
#Reported 148230 pairwise alignments, 148230 HSPs.
#19558 queries aligned.
makehomologs.py -i Tmi_V4b17_hintsutr_vs_TwiV4_AUG.blastp.n.tab -f augustus.hints_utr.aa.nr.n.fasta TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.prot.n.fasta -p 234 -o tethya_clusters_v1 -z 2 -s 2 -M 10 -H 9 -c

for FILE in clusters_tethya_clusters_v1/*.fasta ; do BASE="${FILE%.fasta}" ; mafft $FILE > $BASE.aln ; done

for FILE in clusters_tethya_clusters_v1/homologs_*.aln ; do alignment_conserved_site_to_dots.py -t -n -a $FILE >> Tmi_V4b17_hintsutr_vs_TwiV4_AUG.homologs_identity.tab ; done
```



