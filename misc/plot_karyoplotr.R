#load in data
setwd('~/Desktop/')

genomefile <- read.table("~/project storage/project_comparative_genome/genomes/davidsonii/annot_Pdavidsonii_1mb_chromo_info.txt")
colnames(genomefile) <- c("chr", "start", "end")

annofile <- read.table("~/project storage/project_comparative_genome/genomes/davidsonii/annot_Pdavidsonii_1mb.gffread.genes.bed")
annofile <- annofile[,c(1,2,3,6,4)]
colnames(annofile) <- c("chr", "start", "end", "strand", "gene_id")


#For plotting genic content along scaffolds with karyoploteR:
########################################
library(karyoploteR)

#convert to Granges object
gr <- GenomicRanges::makeGRangesFromDataFrame(annofile)





#plot genome karyotypes and add windowed CDS density
pdf("CDS_density_100kbwindows.pdf")
kp <- plotKaryotype(genome = genomefile, labels.plotter = NULL)
kpAddChromosomeNames(karyoplot = kp, cex = 0.5)

kpPlotDensity(karyoplot = kp,
              data = gr,
              window.size = 100000, col = "black")
dev.off()



#here's a different visualization
pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

pdf("CDS_density_1Mbwindows_view4.pdf")
kp <- plotKaryotype(genome=genomefile, plot.type=4, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp,
                    main="Gene Density")
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45, cex = 0.5)
kpPlotDensity(kp, gr, window.size = 1000000, col="grey")




#plot genome karyotypes and add windowed CDS density
pdf("CDS_density_100kbwindows_view6.pdf")
kp <- plotKaryotype(genome = genomefile, plot.type=6, main = "Gene Density",
                    ideogram.plotter = NULL, cex = 0.5)
kp <- kpPlotDensity(kp, gr, window.size = 100000, data.panel="ideogram", col="orange", border="orange", r0=0.5, r1=1)
kp <- kpPlotDensity(kp, gr, window.size = 100000, data.panel="ideogram", col="orange", border="orange", r0=0.5, r1=0)
dev.off()
########################################




