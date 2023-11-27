# plot T wilhelma contig gc coverage and color by bacterial bin
# by WRF last modified 2023-04-07

inputfile = "~/git/Tethya_wilhelma_genome/02-assembly/Tethya_wilhelma_old_versions/v4_10x/validation/Tethya_wilhelma_V4_P_RNA_scaffold.hisat2_dna.gc_cov.tab"

coveragedata = read.table(inputfile, header=TRUE, sep='\t')
coveragedata = coveragedata[rev(1:nrow(coveragedata)),]

scafnames = coveragedata[,1]

contiglengths = coveragedata[["length"]]
magnituderange = range(log(contiglengths,base=10))
# size of each point
pchsize = log(contiglengths,base=10)-magnituderange[1]
# sizes of three reference points in the legend
legendsizes = c(0.25,0.5,0.75) * (magnituderange[2]-magnituderange[1]) + (magnituderange[1])
legendlabels = round(10^legendsizes)
legendpch = legendsizes - magnituderange[1]

totalsize = sum(contiglengths)

# get bin identity
twi_bact_bins = read.table("~/git/Tethya_wilhelma_genome/03-bacteria/bins/twi_scaf_bin_identity.tab", header = FALSE, sep = "\t") 
twi_sponge_scafs = read.table("~/git/Tethya_wilhelma_genome/02-assembly/Tethya_wilhelma_old_versions/v4_10x/validation/Tethya_wilhelma_V4_P_RNA_scaffold.moleculo_minimap2.b47.names")
bin_id = twi_bact_bins$V2[match(scafnames, twi_bact_bins[["V1"]] )]
bin_id[is.na(bin_id)] = 0
bin_id[match(twi_sponge_scafs[["V1"]], scafnames)] = -1

# assign to colors
#               green       gray unassign    purple      pink          red       red          red          red
pointcolor = c("#10981c33", "#00000044", "#4801939a", "#de68909a", "#b4003dc1", "#b4003dc1", "#b4003dc1", "#b4003dc1" )

# generate figure
outputfile = "~/git/Tethya_wilhelma_genome/03-bacteria/figures/Tethya_wilhelma_V4.metabat_bins.pdf"
pdf(file=outputfile, width=8, height=6)
par(mar=c(4.5,4.5,3,1))
plot(coveragedata[["coverage"]], coveragedata[["GC"]], type='p', 
     xlim=c(0,530), ylim=c(25,75), frame.plot=FALSE, 
     xlab="Mean coverage of mapped DNA reads", ylab="GC%", main="Twi assembly V4 and bacterial bins",
     pch=16, col=pointcolor[(bin_id+2)], cex=pchsize, 
     cex.axis=1.5, cex.lab=1.4 )
legend(400,75, legend=c("Bin 1 (7.58Mb)", "Bin 2 (3.48Mb)", "Bin 3-6"), pch=15, pt.cex=2, bty = 'n',
       col=c("#480193cc", "#de6890cc", "#b4003dcc"), cex=1.1)
legend(400,30, legend=c("Sponge"), pch=15, pt.cex=2, bty = 'n',
       col=c("#10981ccc"), cex=1.1)
dev.off()

png(file="~/git/Tethya_wilhelma_genome/03-bacteria/figures/Tethya_wilhelma_V4.metabat_bins.png", width=8, height=6, units = "in", res=90)
par(mar=c(4.5,4.5,3,1))
plot(coveragedata[["coverage"]], coveragedata[["GC"]], type='p', 
     xlim=c(0,530), ylim=c(25,75), frame.plot=FALSE, 
     xlab="Mean coverage of mapped DNA reads", ylab="GC%", main="Twi assembly V4 and bacterial bins",
     pch=16, col=pointcolor[(bin_id+2)], cex=pchsize, 
     cex.axis=1.5, cex.lab=1.4 )
legend(400,75, legend=c("Bin 1 (7.58Mb)", "Bin 2 (3.48Mb)", "Bin 3-6"), pch=15, pt.cex=2, bty = 'n',
       col=c("#480193cc", "#de6890cc", "#b4003dcc"), cex=1.1)
legend(400,30, legend=c("Sponge"), pch=15, pt.cex=2, bty = 'n',
       col=c("#10981ccc"), cex=1.1)
dev.off()

################################################################################

inputfile = "~/git/Tethya_wilhelma_genome/02-assembly/Tethya_minuta/blobtools/testblob_modif.blobDB.table.txt"

blobplot_headers = c("name", "length", "GC", "N", "bam0", "phylum.t.6", "phylum.s.7", "phylum.c.8")
coveragedata = read.table(inputfile, header=FALSE, sep='\t', col.names = blobplot_headers)
#coveragedata = coveragedata[rev(1:nrow(coveragedata)),]

tmi_bact_bins = read.table("~/git/Tethya_wilhelma_genome/02-assembly/Tethya_minuta/TmiV4_metabat_bin.clusters.tab")
bin_number_id = names(table(tmi_bact_bins$V2))
bin_ids = match(tmi_bact_bins$V2, bin_number_id)
sum(table(tmi_bact_bins$V2)[2:24])
# pointcolor = c("#6890de9a", "#b4003dc1", "#10981c33", "#10981c33", "#b4003dc1", "#b4003dc1", 
#                "#003db4c1", "#10981c33", "#de68909a", "#10981c33", "#10981c33", "#10981c33", 
#                "#6890de9a", "#b4003dc1", "#10981c33", "#10981c33", "#003db4c1", "#10981c33", 
#                "#6890de9a", "#00000044", "#6890de9a", "#10981c33", "#4801939a", "#de68909a" )

tmi_bact_bins
#> names(table(coveragedata$phylum.t.6))
#[1] "Acidobacteria"  "Actinobacteria" "Annelida"       "Arthropoda"     "Bacteria-undef" "Bacteroidetes" 
#[7] "Basidiomycota"  "Brachiopoda"    "Chlamydiae"     "Chordata"       "Cnidaria"       "Echinodermata" 
#[13] "Euryarchaeota"  "Firmicutes"     "Hemichordata"   "Mollusca"       "Mucoromycota"   "Nematoda"      
#[19] "Nitrospirae"    "no-hit"         "Planctomycetes" "Porifera"       "Proteobacteria" "Uroviricota"

#bin_id = match( coveragedata$phylum.t.6, names(table(coveragedata$phylum.t.6)))
bin_id = match( coveragedata$name, tmi_bact_bins$V1 )[bin_ids]
scafnames = coveragedata[,1]
contiglengths = coveragedata[["length"]]
magnituderange = range(log(contiglengths,base=10))
pchsize = log(contiglengths,base=10)-magnituderange[1]

#   1  2   3   4   5  6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
#   0  14  31  45  80 102 105 116 145 149 237 241 260 280 281 302 485 502 560 564 743 758 764 808 
# 239   8   3 255   8   1  12   8   9   6   7   8  54  95   9  40  36  14  40  29  30  20 110   2 

pointcolor = c("#00000044", "green", "green", "#10981c33", "yellow", "purple", 
               "#045a8d88", "black", "orange", "red", "blue", "#de68909a", 
               "pink", "#ce125688", "red", "purple", "green", "red", 
               "darkgreen", "black", "pink", "#12459844", "#679abd", "#fbc961" )

#pdf("~/git/Tethya_wilhelma_genome/02-assembly/Tethya_minuta/tminuta_V4_bins_noID.pdf", width = 8, height = 6)
plot(coveragedata[["bam0"]], coveragedata[["GC"]]*100, type='p', 
     xlim=c(0,50), ylim=c(25,75), frame.plot=FALSE, 
     xlab="Mean coverage of mapped DNA reads", ylab="GC%", main="Tmi assembly V4 and bacterial bins",
     pch=16, col=pointcolor[bin_id], cex=pchsize, 
     cex.axis=1.5, cex.lab=1.4 )
#dev.off()


#