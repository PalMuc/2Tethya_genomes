#!/bin/bash

# download all SRA from Mills 2018
# https://elifesciences.org/articles/31176
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA380886
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395998 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395997 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395986 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395985 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395984 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395983 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395982 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395981 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395980 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395979 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395978 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395977 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395976 
~/sratoolkit.2.11.3-ubuntu64/bin/fastq-dump --gzip --defline-seq '@$sn[_$rn]/$ri' SRR5395975

# rename by treatment factor
mv  SRR5395998.fastq.gz LowO2_Treat_5_SRR5395998.fastq.gz
mv  SRR5395997.fastq.gz LowO2_Treat_4_SRR5395997.fastq.gz
mv  SRR5395986.fastq.gz LowO2_Treat_3_SRR5395986.fastq.gz
mv  SRR5395985.fastq.gz LowO2_Treat_2_SRR5395985.fastq.gz
mv  SRR5395984.fastq.gz LowO2_Treat_1_SRR5395984.fastq.gz
mv  SRR5395983.fastq.gz LowO2_Control_3_SRR5395983.fastq.gz
mv  SRR5395982.fastq.gz LowO2_Control_2_SRR5395982.fastq.gz
mv  SRR5395981.fastq.gz LowO2_Control_1_SRR5395981.fastq.gz
mv  SRR5395980.fastq.gz Draw-down-GW3644_SRR5395980.fastq.gz
mv  SRR5395979.fastq.gz Draw-down-GW3643_SRR5395979.fastq.gz
mv  SRR5395978.fastq.gz Draw-down-GW3642_SRR5395978.fastq.gz
mv  SRR5395977.fastq.gz Anoxia_shock_GW3641-1_SRR5395977.fastq.gz
mv  SRR5395976.fastq.gz Anoxia_shock_GW3639_SRR5395976.fastq.gz
mv  SRR5395975.fastq.gz Anoxia_shock_GW3640_SRR5395975.fastq.gz

# combine and run hisat
{ ~/hisat2-2.2.1/hisat2 -p 6 -q -x Tethya_wilhelma_V4_P_RNA_scaffold.fasta -U Anoxia_shock_GW3639_SRR5395976.fastq.gz,Anoxia_shock_GW3640_SRR5395975.fastq.gz,Anoxia_shock_GW3641-1_SRR5395977.fastq.gz 2>&3 | ~/samtools-1.14/samtools sort - -o Anoxia_shock_combined_v4.bam ; } 3> Anoxia_shock_combined_v4.log
{ ~/hisat2-2.2.1/hisat2 -p 6 -q -x Tethya_wilhelma_V4_P_RNA_scaffold.fasta -U Draw-down-GW3642_SRR5395978.fastq.gz,Draw-down-GW3643_SRR5395979.fastq.gz,Draw-down-GW3644_SRR5395980.fastq.gz 2>&3 | ~/samtools-1.14/samtools sort - -o Draw-down_combined_v4.bam ; } 3> Draw-down_combined_v4.log
{ ~/hisat2-2.2.1/hisat2 -p 6 -q -x Tethya_wilhelma_V4_P_RNA_scaffold.fasta -U LowO2_Control_1_SRR5395981.fastq.gz,LowO2_Control_2_SRR5395982.fastq.gz,LowO2_Control_3_SRR5395983.fastq.gz 2>&3 | ~/samtools-1.14/samtools sort - -o LowO2_Control_combined_v4.bam ; } 3> LowO2_Control_combined_v4.log
{ ~/hisat2-2.2.1/hisat2 -p 6 -q -x Tethya_wilhelma_V4_P_RNA_scaffold.fasta -U LowO2_Treat_2_SRR5395985.fastq.gz,LowO2_Treat_3_SRR5395986.fastq.gz,LowO2_Treat_4_SRR5395997.fastq.gz,LowO2_Treat_5_SRR5395998.fastq.gz 2>&3 | ~/samtools-1.14/samtools sort - -o LowO2_Treat_combined_v4.bam ; } 3> LowO2_Treat_combined_v4.log

# convert bams to beds, and beds to bigwigs
for FILE in *v4.bam ; do BASE="${FILE%.bam}" ; ~/bedtools2/bin/bedtools genomecov -ibam $BASE.bam -bga -split > $BASE.bed ; ~/ucsc_genome_tools/bedSort $BASE.bed $BASE.bed ; ~/ucsc_genome_tools/bedGraphToBigWig $BASE.bed Tethya_wilhelma_V4_P_RNA_scaffold.sizes $BASE.bw ; done

# run salmon on the transcripts
# https://combine-lab.github.io/salmon/faq/
# https://salmon.readthedocs.io/en/latest/salmon.html
cat TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.no_comment.nucl.fasta Tethya_wilhelma_V4_P_RNA_scaffold.rnum.fasta > TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.w_decoys.fasta
~/salmon-latest_linux_x86_64/bin/salmon index -t TwilhelmaV4_AUGUSTUS_predicted_genes_below_50_perc.rnum.w_decoys.fasta -i TwiV4_AUG_tx_index -d decoys.txt -k 31 -p 4

~/salmon-latest_linux_x86_64/bin/salmon quant -i TwiV4_AUG_tx_index -l SR -r Anoxia_shock_GW3639_SRR5395976.fastq.gz Anoxia_shock_GW3640_SRR5395975.fastq.gz Anoxia_shock_GW3641-1_SRR5395977.fastq.gz --validateMappings -o Anoxia_shock_combined_v4.s.quant

~/salmon-latest_linux_x86_64/bin/salmon quant -i TwiV4_AUG_tx_index -l SR -r Draw-down-GW3642_SRR5395978.fastq.gz Draw-down-GW3643_SRR5395979.fastq.gz Draw-down-GW3644_SRR5395980.fastq.gz --validateMappings -o Draw_down_combined_v4.s.quant

~/salmon-latest_linux_x86_64/bin/salmon quant -i TwiV4_AUG_tx_index -l SR -r LowO2_Control_1_SRR5395981.fastq.gz LowO2_Control_2_SRR5395982.fastq.gz LowO2_Control_3_SRR5395983.fastq.gz --validateMappings -o LowO2_Control_combined_v4.s.quant

~/salmon-latest_linux_x86_64/bin/salmon quant -i TwiV4_AUG_tx_index -l SR -r LowO2_Treat_1_SRR5395984.fastq.gz LowO2_Treat_2_SRR5395985.fastq.gz LowO2_Treat_3_SRR5395986.fastq.gz LowO2_Treat_4_SRR5395997.fastq.gz LowO2_Treat_5_SRR5395998.fastq.gz --validateMappings -o LowO2_Treat_combined_v4.s.quant




#
