### This is the data repository for the following publication 

## The genomes of the aquarium sponges *Tethya wilhelma* and *Tethya minuta* (Porifera: Demospongiae)

### Wörheide, Gert <sup>1,2,3,*</sup>, Francis, Warren R.<sup>1</sup>; Deister, Fabian<sup>1</sup>; Stefan Krebs<sup>4</sup>, Erpenbeck, Dirk<sup>1,2</sup>; Vargas, Sergio<sup>1</sup>

<sup>1</sup> Department of Earth and Environmental Sciences, Paleontology and Geobiology, Ludwig-Maximilians-Universität München, Munich, Germany.<br>
<sup>2</sup> GeoBio-Center, Ludwig-Maximilians-Universität München, Munich, Germany<br>
<sup>3</sup> Staatliche Naturwissenschaftliche Sammlungen Bayerns (SNSB)–Bayerische Staatssammlung für Paläontologie und Geologie, Munich, Germany<br>
<sup>4</sup> Laboratory for Functional Genome Analysis (LAFUGA), Gene Center, Ludwig-Maximilians-Universität München, Munich, Germany<br>

*corresponding author

### ABSTRACT <br>
Sponges (phylum Porifera) are aquatic sessile metazoans found worldwide in marine and freshwater environments. They are significant in the animal tree of life as one of the earliest-branching metazoan lineages and as filter feeders play crucial ecological roles, particularly in coral reefs, but are susceptible to the effects of climate change. In the face of the current biodiversity crisis, genomic data is crucial for species conservation efforts and predicting their evolutionary potential in response to environmental changes. However, there is a limited availability of culturable sponge species with annotated high-quality genomes to further comprehensive insights into animal evolution, function, and their response to the ongoing global change. Despite the publication of a few high-quality annotated sponge genomes, there remains a gap in resources for culturable sponge species. To address this gap, we provide high quality draft genomes of the two congeneric aquarium species *Tethya wilhelma* and *Tethya minuta*, small ball-shaped demosponges that are easily maintained long-term in ex situ culture. As such, they offer promising opportunities as laboratory models to contribute to advancing our understanding of sponge biology and provide valuable resources for studying animal evolution, function, and responses to environmental challenges.


### NOTE <br>
*Tethya wilhelma* and *Tethya minuta* are ball-shaped demosponge that undergoes cyclic contractions and are laboratory models for many topics, including multicellularity, early-animal evolution, biomineralization, and microbial interactions. Both were described only in 2001 after being found as bystanders in public aquaria (see [Sara et al 2001](https://www.researchgate.net/publication/215665065_Three_New_Species_of_Tethya_Porifera_Demospongiae_from_German_Aquaria)).

### Please cite:<br>
Wörheide G, Francis WR, Deister F et al. The genomes of the aquarium sponges Tethya wilhelma and Tethya minuta (Porifera: Demospongiae) [version 1; peer review: awaiting peer review].
[ F1000Research 2024, 13:679](http://dx.doi.org/10.12688/f1000research.150836.1)

---

# [01 - preprocessing steps](https://github.com/PalMuc/Tethya_wilhelma_genome/tree/main/01-preprocessing) #

![graphs of albacore readout](https://github.com/PalMuc/Tethya_wilhelma_genome/blob/main/01-preprocessing/Tethya_minuta_read_QC/20180719_1009_Tethya_minuta_FD_1-albacore-2.3.1-FLO-MIN106-SQK-LSK109.fastq-lrplots.png)

# [02 - assembly steps](https://github.com/PalMuc/Tethya_wilhelma_genome/tree/main/02-assembly) #
All bioinformatic steps and versions of assemblies

![barplot of Tethya genome assemblies](https://github.com/PalMuc/Tethya_wilhelma_genome/raw/main/02-assembly/figures/Tethya_sp_combined.sizes.png)

# [03 - bacterial contig separation](https://github.com/PalMuc/Tethya_wilhelma_genome/tree/main/03-bacteria) #
Code for identification and separation of contigs from associated bacteria, for both *Tethya wilhelma* and *Tethya minuta*

![blob plot of Tethya wilhelma](https://github.com/PalMuc/Tethya_wilhelma_genome/raw/main/03-bacteria/figures/Tethya_wilhelma_V4.metabat_bins.png)

# [04 - mitochondrial genomes](https://github.com/PalMuc/Tethya_wilhelma_genome/tree/main/04-mitochondria) #
All 4 complete genomes are rotated to begin with COX1. All genes are on the same (forward) strand, as found in many other sponges.

![schematic of Tethya mito genomes](https://github.com/PalMuc/Tethya_wilhelma_genome/raw/main/04-mitochondria/figures/MITOS_output_graph.png)

# [05 - gene annotation](https://github.com/PalMuc/Tethya_wilhelma_genome/tree/main/05-annotation) #
Annotation files for V4 of *Tethya wilhelma* and V4 of *Tethya minuta*, both using AUGUSTUS.


---

Do **NOT** cite the unpublished [2017 bioRxiv preprint](https://www.biorxiv.org/content/10.1101/120998v3), which was `v1`.

This project was seed-funded by the LMUexcellent program (Project MODELSPONGE) to G.W. and D.E. through the German Excellence Initiative, and also benefited from funding by VILLUM FONDEN (Grant 9278 "Early evolution of multicellular sponges") to G.W.

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.



