## DATA: T-citrina_nanopore_assembly.fasta

```shell
longest scaffold = 644682
shortest scaffold = 771
N50 = 59128
number of sequences = 1593
number of N|n-position = 1927330
N|n percent = 3.36453194368497
total length = 57283748
```

```shell
#!/bin/bash
#
#SBATCH --job-name=augustus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000

./../program/augustus-3.3.1/bin/augustus --strand=both --codingseq=on --genemodel=complete \
--AUGUSTUS_CONFIG_PATH=/naslx/projects/pn69xe/di52zuc/program/augustus-3.3.1/config \
--species=Tethya_wilh --outfile=T-citrina_predicted_genes_coding_sequences.gff \
T-citrina_nanopore_assembly.fasta
```
# gff_export_cds_and_prot.pl :
```shell
13813 predicted proteins
coding_regions_T-citrina_predicted_genes_coding_sequences.gff.fasta
prot_T-citrina_predicted_genes_coding_sequences.gff.fasta
```