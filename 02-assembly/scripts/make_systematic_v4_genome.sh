#!/bin/bash

################################################################################
################################################################################
# using the manually selected scaffolds from the interactive viewer
# https://github.com/wrf/lavaLampPlot

cut -f 1 Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.tab > Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names
getAinB.py Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names Tethya_wilhelma_V4_P_RNA_scaffold.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.b47.fasta
excludeAinB.py Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names Tethya_wilhelma_V4_P_RNA_scaffold.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.a47.fasta

# rename in size order after scaffolding
~/git/genomeGTFtools/number_contigs_by_length.py Tethya_wilhelma_V4_P_RNA_scaffold.b47.fasta -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -n TwiV4 -o -d . > Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta
gzip Tethya_wilhelma_V4_P_RNA_scaffold.fasta
sizecutter.py -f Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta > Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes

################################################################################
################################################################################
# change all GFFs

# rename scaffold names in AUGUSTUS gene models
~/git/genomeGTFtools/rename_gtf_contigs.py --omit-comments -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g predicted_genes_Tethya_wilhelma_V4_below_50_perc.gff.gz > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff
# Converted 406745 lines and could not change 158

#extract_augustus_features.py -n TwiV4 -c TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.nucl.fasta -p TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.prot.fasta TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff
# Extracting features from TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff  Wed Jan 25 18:02:37 2023
# Finished parsing  Wed Jan 25 18:02:39 2023
# No CDS found in TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff
# Counted 28113 proteins

~/git/genomeGTFtools/misc/augustus_to_gff3.py TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.gff > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.gff


################################################################################
################################################################################
# change all BED files to new sysematic numbering

# rename all bed files, and convert to bigwigs
~/git/genomeGTFtools/rename_gtf_contigs.py -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g Anoxia_shock_combined_v4.bed > Anoxia_shock_combined_v4.rnum.bed
# Converted 12206198 lines and could not change 12079

~/git/genomeGTFtools/rename_gtf_contigs.py -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g Draw-down_combined_v4.bed > Draw_down_combined_v4.rnum.bed
# Converted 12424658 lines and could not change 17032

~/git/genomeGTFtools/rename_gtf_contigs.py -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g LowO2_Control_combined_v4.bed > LowO2_Control_combined_v4.rnum.bed
# Converted 13788021 lines and could not change 119814

~/git/genomeGTFtools/rename_gtf_contigs.py -n -c Tethya_wilhelma_V4_P_RNA_scaffold.rnum.vector -g LowO2_Treat_combined_v4.bed > LowO2_Treat_combined_v4.rnum.bed
# Converted 15223818 lines and could not change 44773

# zip old bg and remove old bw files
gzip Anoxia_shock_combined_v4.bed ; gzip Draw-down_combined_v4.bed ; gzip LowO2_Control_combined_v4.bed ; gzip LowO2_Treat_combined_v4.bed
rm *v4.bw

# 
~/ucsc_genome_tools/bedSort Anoxia_shock_combined_v4.rnum.bed Anoxia_shock_combined_v4.rnum.bed
~/ucsc_genome_tools/bedSort Draw_down_combined_v4.rnum.bed Draw_down_combined_v4.rnum.bed
~/ucsc_genome_tools/bedSort LowO2_Control_combined_v4.rnum.bed LowO2_Control_combined_v4.rnum.bed
~/ucsc_genome_tools/bedSort LowO2_Treat_combined_v4.rnum.bed LowO2_Treat_combined_v4.rnum.bed

~/ucsc_genome_tools/bedGraphToBigWig Anoxia_shock_combined_v4.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes Anoxia_shock_combined_v4.rnum.bw
~/ucsc_genome_tools/bedGraphToBigWig Draw_down_combined_v4.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes Draw_down_combined_v4.rnum.bw
~/ucsc_genome_tools/bedGraphToBigWig LowO2_Control_combined_v4.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes LowO2_Control_combined_v4.rnum.bw
~/ucsc_genome_tools/bedGraphToBigWig LowO2_Treat_combined_v4.rnum.bed Tethya_wilhelma_V4_P_RNA_scaffold.rnum.sizes LowO2_Treat_combined_v4.rnum.bw

for FILE in *.rnum.bed ; do gzip $FILE ; done

################################################################################
################################################################################



