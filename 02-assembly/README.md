# 02 - assemblies #

![barplot comparing assemblies](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/02-assembly/figures/Tethya_sp_combined.sizes.png)


## assembly gaps ##

```
zcat v1_moleculo/tethya-0_1.fa.gz | ~/git/genomeGTFtools/repeat2gtf.py -l - | gzip > v1_moleculo/tethya-0_1.n_gaps.gff
# Parsing repeats of N and n from <stdin>   12:22:36 2023
# Counted 7947 sequences   12:22:44 2023
# Counted 4300 repeats of 1539645 total bases, average 358.06
# Longest repeat was 6774 bases on tethya01_919

zcat v2_dovetail/sponge_tethya_wilhelma_26Apr2018_OM745.fasta.gz | ~/git/genomeGTFtools/repeat2gtf.py -l - | gzip > v2_dovetail/sponge_tethya_wilhelma_26Apr2018_OM745.n_gaps.gff.gz
# Parsing repeats of N and n from <stdin>   12:24:15 2023
# Counted 1353 sequences   12:24:23 2023
# Counted 11248 repeats of 2252241 total bases, average 200.23
# Longest repeat was 6774 bases on Scaffold_1346;HRSCAF=1752

zcat v3_nanopore/gaps_closed_scaffolded_10_x_raw.fasta.gz | ~/git/genomeGTFtools/repeat2gtf.py -l - | gzip > v3_nanopore/gaps_closed_scaffolded_10_x_raw.n_gaps.gff.gz
# Parsing repeats of N and n from <stdin>   12:26:15 2023
# Counted 967 sequences   12:26:24 2023
# Counted 2416 repeats of 1069943 total bases, average 442.86
# Longest repeat was 8380 bases on scaffold411|size26987

zcat v4_10x/Tethya_wilhelma_V4_P_RNA_scaffold.fasta.gz | ~/git/genomeGTFtools/repeat2gtf.py -l - | gzip > v4_10x/Tethya_wilhelma_V4_P_RNA_scaffold.n_gaps.gff.gz
# Parsing repeats of N and n from <stdin>   12:28:38 2023
# Counted 891 sequences   12:28:48 2023
# Counted 2492 repeats of 1077543 total bases, average 432.40
# Longest repeat was 8380 bases on scaffold411


```

## making working versions of AUGUSTUS annotations
Using accessory script [augustus_to_gff3.py](https://github.com/wrf/genomeGTFtools/blob/master/misc/augustus_to_gff3.py)

```
~/git/genomeGTFtools/misc/augustus_to_gff3.py TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff

~/git/genomeGTFtools/misc/augustus_to_gff3.py Tethya_minuta_V4.gff > Tethya_minuta_V4.no_comment.gff

~/git/genomeGTFtools/misc/augustus_to_gff3.py T_citrina_nanopore_assembly_AUGUSTUS.gff > T_citrina_nanopore_assembly_AUGUSTUS.no_comment.gff
```



