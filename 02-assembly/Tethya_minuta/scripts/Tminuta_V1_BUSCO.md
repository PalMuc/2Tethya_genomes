
## Data: *Tethya minuta* assembly WTDBG2. Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa

```shell
#!/bin/bash
#
#SBATCH --job-name=busco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=f.deister@palmuc.org
#SBATCH --mail-type=ALL
#SBATCH --clusters=mpp2
#SBATCH --time=5:00:00
#SBATCH --mem=55000


#export AUGUSTUS_CONFI_PATH=/naslx/projects/pn69xe/di52zuc/program/augustus-3.3.1/config
export PATH=$PATH:/naslx/projects/pn69xe/di52zuc/program/augustus-3.3.1/config
export PATH=$PATH:/naslx/projects/pn69xe/di52zuc/program/augustus-3.3.1/bin
 
python ../scripts/run_BUSCO.py --in Tminuta_assembly_WTDBG_nanopore.ctg.3rd.fa -o OUTPUT -l /naslx/projects/pn69xe/di52zuc/program/busco-master/metazoa_odb9 \
-sp amphimedon  -c 28 -f -m geno >logfile 2>error
```

```shell
INFO    Training Augustus using Single-Copy Complete BUSCOs:
INFO    Converting predicted genes to short genbank files at 04/26/2019 13:18:42...
INFO    All files converted to short genbank files, now running the training scripts at 04/26/2019 13:18:55...
INFO    Pre-Augustus scaffold extraction...
INFO    Re-running Augustus with the new metaparameters, number of target BUSCOs: 474
INFO    [augustus]      47 of 468 task(s) completed at 04/26/2019 13:19:01
INFO    [augustus]      94 of 468 task(s) completed at 04/26/2019 13:19:05
INFO    [augustus]      188 of 468 task(s) completed at 04/26/2019 13:19:12
INFO    [augustus]      234 of 468 task(s) completed at 04/26/2019 13:19:16
INFO    [augustus]      281 of 468 task(s) completed at 04/26/2019 13:19:19
INFO    [augustus]      328 of 468 task(s) completed at 04/26/2019 13:19:29
INFO    [augustus]      375 of 468 task(s) completed at 04/26/2019 13:19:35
INFO    [augustus]      422 of 468 task(s) completed at 04/26/2019 13:19:39
INFO    [augustus]      468 of 468 task(s) completed at 04/26/2019 13:19:48
INFO    Extracting predicted proteins...
INFO    ****** Step 3/3, current time: 04/26/2019 13:19:55 ******
INFO    Running HMMER to confirm orthology of predicted proteins:
INFO    [hmmsearch]     464 of 464 task(s) completed at 04/26/2019 13:19:55
INFO    Results:
INFO    C:65.7%[S:63.2%,D:2.5%],F:10.0%,M:24.3%,n:978
INFO    642 Complete BUSCOs (C)
INFO    618 Complete and single-copy BUSCOs (S)
INFO    24 Complete and duplicated BUSCOs (D)
INFO    98 Fragmented BUSCOs (F)
INFO    238 Missing BUSCOs (M)
INFO    978 Total BUSCO groups searched
INFO    BUSCO analysis done with WARNING(s). Total running time: 327.6717665195465 seconds
INFO    Results written in /naslx/projects/pn69xe/di52zuc/program/busco-master/minuta_Project2/run_OUTPUT/
```