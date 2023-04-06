*De Novo* assembly of *Tethya citrina* using nanopore reads
```shell
#!/bin/bash
#
#SBATCH --job-name=WTDBG
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=48:00:00
#SBATCH --mem=55000


wtdbg2 -x ont -g 100000000 -i ../20190315_T_citrina_PCR-guppy-2.3.5.porechop.fastq -t112 -fo T_citrina_PCR-guppy-2.3.5

wtpoa-cns -t112 -i T_citrina_PCR-guppy-2.3.5.ctg.lay.gz -fo T_citrina_PCR-guppy-2.3.5.ctg.fa

minimap2 -t 112 -x map-pb \
-a T_citrina_PCR-guppy-2.3.5.ctg.fa ../20190315_T_citrina_PCR-guppy-2.3.5.porechop.fastq | samtools view -Sb  >T_citrina_PCR-guppy-2.3.5.ctg.map.bam

samtools sort T_citrina_PCR-guppy-2.3.5.ctg.map.bam > T_citrina_PCR-guppy-2.3.5.ctg.map.srt.bam
samtools view T_citrina_PCR-guppy-2.3.5.ctg.map.srt.bam | wtpoa-cns -t 112 -d T_citrina_PCR-guppy-2.3.5.ctg.fa \
-i - -fo T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa

bwa mem -t 112 T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa \
link_Galaxy23-Trimmomatic_on_Taurantium_Ind01_Lib02_gDNA_Read1.fastq_R1_paired.fastq \
link_Galaxy24-Trimmomatic_on_Taurantium_Ind01_Lib02_gDNA_Read2.fastq_R2_paired.fastq \
| samtools view -Sb - >sr.bam
samtools sort sr.bam >sr.srt.bam
samtools view sr.srt.bam | wtpoa-cns -t 16 -x sam-sr -d T_citrina_PCR-guppy-2.3.5.ctg.2nd.fa -i - -fo T_citrina_PCR-guppy-2.3.5.ctg.3rd.fa


```

```shell
--
-- total memory       65256936.0 kB
-- available          63126980.0 kB
-- 28 cores
-- Starting program: wtdbg2 -x ont -p 0 -S 2 -L 1000 -s 0.1 -l 1000 -m 100 -g 110m -i T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.fastq -fo T_citrina_27_03_2019_FD-guppy-2.3.7.porechop
-- pid                      8133
-- date         Fri Apr 12 17:02:32 2019
--
[Fri Apr 12 17:02:32 2019] loading reads
534963 reads
[Fri Apr 12 17:03:02 2019] Done, 534963 reads (>=1000 bp), 1489594588 bp, 5557763 bins
** PROC_STAT(0) **: real 30.421 sec, user 6.560 sec, sys 3.160 sec, maxrss 524720.0 kB, maxvsize 714236.0 kB
[Fri Apr 12 17:03:02 2019] Set --edge-cov to 2
KEY PARAMETERS: -k 15 -p 0 -K 1000.049988 -A -S 2.000000 -s 0.100000 -g 110000000 -X 50.000000 -e 2 -L 1000
[Fri Apr 12 17:03:02 2019] generating nodes, 4 threads
[Fri Apr 12 17:03:02 2019] indexing bins[0,5557763] (1422787328 bp), 4 threads
[Fri Apr 12 17:03:02 2019] - scanning kmers (K15P0S2.00) from 5557763 bins
5557763 bins
********************** Kmer Frequency **********************
                                                                                                    
                                                                                                    
|                                                                                                   
|                                                                                                   
||                                                                                                  
||                                                                                                  
||                                                                                                  
||                                                                                                  
||                                                                                                  
|||                                                                                                 
|||                                                                                                 
||||                                                                                                
||||                                                                                                
|||||                                                                                               
|||||                                                                                               
||||||                                                                                              
|||||||                                                                                             
||||||||                                                                                            
||||||||||                                                                                          
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
**********************     1 - 201    **********************
Quatiles:
   10%   20%   30%   40%   50%   60%   70%   80%   90%   95%
     2     3     4     5     7    10    14    22    58   203
# If the kmer distribution is not good, please kill me and adjust -k, -p, and -K
# Cannot get a good distribution anyway, should adjust -S -s, also -A -e in assembly
** PROC_STAT(0) **: real 79.007 sec, user 164.550 sec, sys 13.170 sec, maxrss 3046280.0 kB, maxvsize 3622948.0 kB
[Fri Apr 12 17:03:51 2019] - high frequency kmer depth is set to 1000
[Fri Apr 12 17:03:51 2019] - Total kmers = 163732398
[Fri Apr 12 17:03:51 2019] - average kmer depth = 5
[Fri Apr 12 17:03:51 2019] - 62883857 low frequency kmers (<2)
[Fri Apr 12 17:03:51 2019] - 4522 high frequency kmers (>1000)
[Fri Apr 12 17:03:51 2019] - indexing 100844019 kmers, 589113020 instances (at most)
5557763 bins
[Fri Apr 12 17:05:07 2019] - indexed  100844019 kmers, 588349178 instances
[Fri Apr 12 17:05:07 2019] - masked 219 bins as closed
[Fri Apr 12 17:05:07 2019] - sorting
** PROC_STAT(0) **: real 158.051 sec, user 448.480 sec, sys 30.500 sec, maxrss 6036220.0 kB, maxvsize 8538244.0 kB
[Fri Apr 12 17:05:10 2019] Done
534963 reads|total hits 22260702
** PROC_STAT(0) **: real 649.679 sec, user 2312.800 sec, sys 106.190 sec, maxrss 6765772.0 kB, maxvsize 9735852.0 kB
[Fri Apr 12 17:13:22 2019] clipping ... 12.96% bases
[Fri Apr 12 17:13:23 2019] generating regs ... 23056250
[Fri Apr 12 17:13:29 2019] sorting regs ...  Done
[Fri Apr 12 17:13:30 2019] generating intervals ...  939464 intervals
[Fri Apr 12 17:13:30 2019] selecting important intervals from 939464 intervals
[Fri Apr 12 17:13:32 2019] Intervals: kept 102450, discarded 837014
** PROC_STAT(0) **: real 660.497 sec, user 2330.960 sec, sys 108.440 sec, maxrss 6765772.0 kB, maxvsize 9735852.0 kB
[Fri Apr 12 17:13:32 2019] Done, 102450 nodes
[Fri Apr 12 17:13:32 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.1.nodes". Done.
[Fri Apr 12 17:13:33 2019] masked 108 high coverage nodes (>200 or <2)
[Fri Apr 12 17:13:33 2019] masked 7308 repeat-like nodes by local subgraph analysis
[Fri Apr 12 17:13:33 2019] generating edges
[Fri Apr 12 17:13:33 2019] Done, 399993 edges
[Fri Apr 12 17:13:33 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.1.reads". Done.
[Fri Apr 12 17:13:34 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.1.dot.gz". Done.
[Fri Apr 12 17:13:34 2019] graph clean
[Fri Apr 12 17:13:34 2019] rescued 3696 low cov edges
[Fri Apr 12 17:13:34 2019] deleted 216 binary edges
[Fri Apr 12 17:13:34 2019] deleted 25117 isolated nodes
[Fri Apr 12 17:13:34 2019] cut 6937 transitive edges
[Fri Apr 12 17:13:34 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.2.dot.gz". Done.
[Fri Apr 12 17:13:35 2019] 822 bubbles; 3159 tips; 29 yarns;
[Fri Apr 12 17:13:35 2019] deleted 440 isolated nodes
[Fri Apr 12 17:13:35 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.3.dot.gz". Done.
[Fri Apr 12 17:13:35 2019] cut 173 branching nodes
[Fri Apr 12 17:13:35 2019] deleted 29 isolated nodes
[Fri Apr 12 17:13:35 2019] building unitigs
[Fri Apr 12 17:13:35 2019] TOT 81814528, CNT 6325, AVG 12936, MAX 529664, N50 25344, L50 925, N90 5376, L90 3559, Min 1024
[Fri Apr 12 17:13:35 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.frg.nodes". Done.
[Fri Apr 12 17:13:35 2019] generating links
[Fri Apr 12 17:13:35 2019] generated 2769 links
[Fri Apr 12 17:13:35 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.frg.dot.gz". Done.
[Fri Apr 12 17:13:35 2019] rescue 0 weak links
[Fri Apr 12 17:13:35 2019] deleted 30 binary links
[Fri Apr 12 17:13:35 2019] cut 33 transitive links
[Fri Apr 12 17:13:35 2019] remove 278 boomerangs
[Fri Apr 12 17:13:35 2019] detached 12 repeat-associated paths
[Fri Apr 12 17:13:35 2019] remove 445 weak branches
[Fri Apr 12 17:13:35 2019] cut 11 tips
[Fri Apr 12 17:13:35 2019] pop 0 bubbles
[Fri Apr 12 17:13:35 2019] cut 0 tips
[Fri Apr 12 17:13:35 2019] output "T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.dot.gz". Done.
[Fri Apr 12 17:13:35 2019] building contigs
[Fri Apr 12 17:13:35 2019] searched 4647 contigs
[Fri Apr 12 17:13:35 2019] Estimated: TOT 79388416, CNT 2898, AVG 27395, MAX 702208, N50 42240, L50 548, N90 12032, L90 1905, Min 5120
[Fri Apr 12 17:13:45 2019] output 2898 contigs
[Fri Apr 12 17:13:45 2019] Program Done
** PROC_STAT(TOTAL) **: real 673.280 sec, user 2372.770 sec, sys 109.430 sec, maxrss 6765772.0 kB, maxvsize 9735852.0 kB
---
--
-- total memory       65256936.0 kB
-- available          63142484.0 kB
-- 28 cores
-- Starting program: wtpoa-cns -t112 -i T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.lay.gz -fo T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.fa
-- pid                      9120
-- date         Fri Apr 12 17:13:45 2019
--
2898 contigs 55859 edges 76793372 bases
** PROC_STAT(TOTAL) **: real 365.363 sec, user 5454.010 sec, sys 120.790 sec, maxrss 7501564.0 kB, maxvsize 14784688.0 kB
---
[M::mm_idx_gen::1.680*1.28] collected minimizers
[M::mm_idx_gen::1.974*1.72] sorted minimizers
[M::main::1.975*1.72] loaded/built the index for 2898 target sequence(s)
[M::mm_mapopt_update::2.107*1.67] mid_occ = 24
[M::mm_idx_stat] kmer size: 19; skip: 10; is_hpc: 1; #seq: 2898
[M::mm_idx_stat::2.226*1.64] distinct minimizers: 8448920 (90.87% are singletons); average occurrences: 1.191; average spacing: 7.642
[M::worker_pipeline::118.302*11.40] mapped 388260 sequences
[M::worker_pipeline::162.235*12.49] mapped 388262 sequences
[M::worker_pipeline::204.298*12.90] mapped 387456 sequences
[M::worker_pipeline::222.848*11.84] mapped 351157 sequences
[M::main] Version: 2.13-r850
[M::main] CMD: minimap2 -t 112 -x map-pb -a T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.fa T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.fastq
[M::main] Real time: 223.056 sec; CPU: 2638.205 sec; Peak RSS: 11.538 GB
[bam_sort_core] merging from 5 files and 1 in-memory blocks...
--
-- total memory       65256936.0 kB
-- available          63113072.0 kB
-- 28 cores
-- Starting program: wtpoa-cns -t112 -d T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.fa -i - -fo T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.2nd.fa
-- pid                     10939
-- date         Fri Apr 12 17:26:19 2019
--
2898 contigs 52067 edges 55121056 bases
** PROC_STAT(TOTAL) **: real 2299.195 sec, user 42685.090 sec, sys 492.430 sec, maxrss 5115664.0 kB, maxvsize 14660136.0 kB
```

```shell
sequence_counter.pl T_citrina_PCR-guppy-2.3.5.ctg.3rd.fa 
Longest scaffold = 864048
Shortest scaffold = 143
N50 = 27600
Number of sequences = 3566
Number of N|n-position = 0
N|n percent = 0
Total sequence length = 64090650

```
## SSPACE-long 
```bash
#!/bin/bash
#
#SBATCH --job-name=n_sspace_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=6:00:00
#SBATCH --mem=55000

export PATH=$PATH:/naslx/projects/pn69xe/di52zuc/program/SSPACE-LongRead_v1-1

perl /naslx/projects/pn69xe/di52zuc/program/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl \
-c T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.2nd.fa \
-b Tcitrina_nanopore_reads_sspace \
-i 80 -t 80 -o 5000 -l 1 -k 1 \
-p long_reads_combined.fa 2> nanopore_reads_errorlog.txt >nanopore_reads_contigs_logfile.txt

```


## SSPACE-STANDARD

```
#!/bin/bash
#
#SBATCH --job-name=n_sspace_long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=6:00:00
#SBATCH --mem=55000

export PATH=$PATH:/naslx/projects/pn69xe/di52zuc/program/SSPACE-LongRead_v1-1



perl /naslx/projects/pn69xe/di52zuc/program/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl \
-c T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.ctg.2nd.fa \
-b Tcitrina_nanopore_reads_sspace \
-i 80 -t 80 -o 5000 -l 1 -k 1 \
-p long_reads_combined.fa 2> nanopore_reads_errorlog.txt >nanopore_reads_contigs_logfile.txt
```
```
perl sequence_counter.pl T-citrina_scaffolded_final.fasta 

                Longest scaffold = 644682
                Shortest scaffold = 771
                N50 = 59128
                Number of sequences = 1593
                Number of N|n-position = 1927330
                N|n percent = 3.36453194368497
                Total sequence length = 57283748

```
