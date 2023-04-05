# analyses of Tethya assembly, and quality control
# by WRF 2023-01-26

setwd("~/git/Tethya_wilhelma_genome/02-assembly/")

################################################################################
# T wilhelma scaffold overview barplot
tw_scaf_size_data = read.table("stats/TwiV4_scaffold.rnum.sizes", sep="\t", row.names = 1)

running_total_bases = cumsum(tw_scaf_size_data[,1])
n50_cutoff = sum(tw_scaf_size_data[["V2"]]) * 0.5
n50_index = min(which(running_total_bases > n50_cutoff))
n90_cutoff = sum(tw_scaf_size_data[["V2"]]) * 0.9
n90_index = min(which(running_total_bases > n90_cutoff))

pdf(file="~/git/Tethya_wilhelma_genome/02-assembly/figures/TwiV4_scaffold.rnum.sizes.pdf", width=8, height=5)
par(mar=c(4.5,4.5,1,1))
b = barplot(tw_scaf_size_data$V2[1:40]/1e6, col="#10981cdd", border = NA, xlab="TwiV4 scaffold ID",
        ylab="Scaffold size (Mbp)", cex.lab=1.3, axes=FALSE,
        names.arg = gsub("TwiV4.","",rownames(tw_scaf_size_data))[1:40],
        las=3, cex.names = 0.7 )
axis(2, at=seq(0,25,5), cex.axis=1.3)
text(50,25,paste("Showing 40 out of", nrow(tw_scaf_size_data), "scaffolds"), pos=2 )
text(50,22,paste("N50:", round(tw_scaf_size_data[n50_index,1]/1e6,digits = 1), "Mbp"), pos=2 )
segments( (n50_index)*1.2+0.1 , 0,
          (n50_index)*1.2+0.1 , max(tw_scaf_size_data[,1]), col = "#00000066", lty = 2)
segments(0,tw_scaf_size_data[n50_index,1]/1e6,
         48,tw_scaf_size_data[n50_index,1]/1e6, col = "#00000066", lty = 2)
text( n50_index*1.2+0.1,tw_scaf_size_data[n50_index,1]/1e6*1.2, "N50", pos = 4)
dev.off()


################################################################################

# T minuta scaffold overview barplot
tm_scaf_size_data = read.table("stats/Tethya_minuta_V4.sizes", sep="\t", row.names = 1)
tm_sorted_sizes = sort(tm_scaf_size_data$V2, decreasing = TRUE, index.return = TRUE)

pdf(file="~/git/Tethya_wilhelma_genome/02-assembly/figures/Tethya_minuta_V4.sizes.pdf", width=8, height=5)
par(mar=c(4.5,4.5,1,1))
b = barplot(tm_sorted_sizes$x[1:60]/1e6, col="#20bb2cdd", border = NA, xlab="Tminuta V4 scaffold ID",
            ylab="Scaffold size (Mbp)", cex.lab=1.3, axes=FALSE,
            names.arg = gsub("Scaffolds_","",rownames(tm_scaf_size_data)[tm_sorted_sizes$ix])[1:60],
            las=3, cex.names = 0.7 )
axis(2, at=seq(0,4,0.5), cex.axis=1.3)
text(70,4,paste("Showing 60 out of", nrow(tm_scaf_size_data), "scaffolds"), pos=2 )
dev.off()

################################################################################

# T citrina scaffold overview barplot
tc_scaf_size_data = read.table("stats/Tethya_citrina_nanopore_assembly.sizes", sep="\t", row.names = 1)
tc_sorted_sizes = sort(tc_scaf_size_data$V2, decreasing = TRUE, index.return = TRUE)

pdf(file="~/git/Tethya_wilhelma_genome/02-assembly/figures/Tethya_citrina_nanopore_assembly.sizes.pdf", width=8, height=5)
par(mar=c(4.5,4.5,1,1))
b = barplot(tc_sorted_sizes$x[1:60]/1e3, col="#91d60cdd", border = NA, xlab="Tcitrina V3 scaffold ID",
            ylab="Scaffold size (kbp)", cex.lab=1.3, axes=FALSE,
            names.arg = seq(1,60,1),
            las=3, cex.names = 0.7 )
axis(2, at=seq(0,600,100), cex.axis=1.3)
text(70,600,paste("Showing 60 out of", nrow(tc_scaf_size_data), "scaffolds"), pos=2 )
dev.off()

################################################################################

# generate combined figure of 3 genomes

pdf(file="~/git/Tethya_wilhelma_genome/02-assembly/figures/Tethya_sp_combined.sizes.pdf", width=8, height=11, paper = "a4")
par(mar=c(4.5,4.5,1,1), mfrow=c(3,1))
b = barplot(tw_scaf_size_data$V2[1:40]/1e6, col="#10981cdd", border = NA, xlab="TwiV4 scaffold ID",
            ylab="Scaffold size (Mbp)", cex.lab=1.3, axes=FALSE,
            names.arg = gsub("TwiV4.","",rownames(tw_scaf_size_data))[1:40],
            las=3, cex.names = 0.7 )
axis(2, at=seq(0,25,5), cex.axis=1.3)
text(50,25,paste("Showing 40 out of", nrow(tw_scaf_size_data), "scaffolds"), pos=2 )
text(50,22,paste("N50:", round(tw_scaf_size_data[n50_index,1]/1e6,digits = 1), "Mbp"), pos=2 )
segments( (n50_index)*1.2+0.1 , 0,
          (n50_index)*1.2+0.1 , max(tw_scaf_size_data[,1]), col = "#00000066", lty = 2)
segments(0,tw_scaf_size_data[n50_index,1]/1e6,
         48,tw_scaf_size_data[n50_index,1]/1e6, col = "#00000066", lty = 2)
text( n50_index*1.2+0.1,tw_scaf_size_data[n50_index,1]/1e6*1.2, "N50", pos = 4)

b = barplot(tm_sorted_sizes$x[1:60]/1e6, col="#20bb2cdd", border = NA, xlab="Tminuta V4 scaffold ID",
            ylab="Scaffold size (Mbp)", cex.lab=1.3, axes=FALSE,
            names.arg = gsub("Scaffolds_","",rownames(tm_scaf_size_data)[tm_sorted_sizes$ix])[1:60],
            las=3, cex.names = 0.7 )
axis(2, at=seq(0,4,0.5), cex.axis=1.3)
text(75,4,paste("Showing 60 out of", nrow(tm_scaf_size_data), "scaffolds"), pos=2 )

b = barplot(tc_sorted_sizes$x[1:60]/1e3, col="#91d60cdd", border = NA, xlab="Tcitrina V3 scaffold ID",
            ylab="Scaffold size (kbp)", cex.lab=1.3, axes=FALSE,
            names.arg = seq(1,60,1),
            las=3, cex.names = 0.7 )
axis(2, at=seq(0,600,100), cex.axis=1.3)
text(75,600,paste("Showing 60 out of", nrow(tc_scaf_size_data), "scaffolds"), pos=2 )
dev.off()









#