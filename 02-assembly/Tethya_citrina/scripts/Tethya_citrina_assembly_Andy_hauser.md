## Data: T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.fastq-wtdbg2-xont-g54M.racon1.fasta

```shell
perl sequence_counter.pl T_citrina_27_03_2019_FD-guppy-2.3.7.porechop.fastq-wtdbg2-xont-g54M.racon1.fasta
longest scaffold = 368812
shortest scaffold = 2800
N50 = 27126
number of sequences = 3399
number of N|n-position = 0
N|n percent = 0
total length = 54678763

```

## Hisat2 mapping Illumina paired end
```shell
# This file was produced by samtools stats (1.8+htslib-1.8) and can be plotted using plot-bamstats
# This file contains statistics for all reads.
# The command line was:  stats scaffolds.fasta.sam
# CHK, Checksum [2]Read Names   [3]Sequences    [4]Qualities
# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)
CHK     b92394f4        e3a8dccd        f80f9142
# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.
SN      raw total sequences:    12158394
SN      filtered sequences:     0
SN      sequences:      12158394
SN      is sorted:      0
SN      1st fragments:  6079197
SN      last fragments: 6079197
SN      reads mapped:   12158394
SN      reads mapped and paired:        12158394        # paired-end technology bit set + both mates mapped
SN      reads unmapped: 0
SN      reads properly paired:  11542730        # proper-pair bit set
SN      reads paired:   12158394        # paired-end technology bit set
SN      reads duplicated:       0       # PCR or optical duplicate bit set
SN      reads MQ0:      528846  # mapped and MQ=0
SN      reads QC failed:        0
SN      non-primary alignments: 1503378
SN      total length:   1054461809      # ignores clipping
SN      bases mapped:   1054461809      # ignores clipping
SN      bases mapped (cigar):   1022086046      # more accurate
SN      bases trimmed:  0
SN      bases duplicated:       0
SN      mismatches:     9271175 # from NM fields
SN      error rate:     9.070836e-03    # mismatches / bases mapped (cigar)
SN      average length: 86
SN      maximum length: 95
SN      average quality:        36.7
SN      insert size average:    562.2
SN      insert size standard deviation: 1591.6
SN      inward oriented pairs:  5649933
SN      outward oriented pairs: 124405
SN      pairs with other orientation:   42770
SN      pairs on different chromosomes: 140331
```
