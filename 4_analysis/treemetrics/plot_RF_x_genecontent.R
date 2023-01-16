#load in data
setwd('~/Desktop/')
library(tidyverse)

genomefile <- read.table("~/project storage/project_comparative_genome/genomes/davidsonii/annot_Pdavidsonii_1mb_chromo_info.txt")
colnames(genomefile) <- c("chr", "start", "end")

annofile <- read.table("~/project storage/project_comparative_genome/genomes/davidsonii/annot_Pdavidsonii_1mb.gffread.genes.bed")
annofile <- annofile[,c(1,2,3,6,4)]
colnames(annofile) <- c("chr", "start", "end", "strand", "gene_id")


##setup for RF estimation with 10kb trees:
##################################################
library(ape)
library(phangorn)
library(tidyverse)
setwd('~/Desktop/')

#read in the species tree and reroot to correct outgroup
speciestree <- read.tree("~/project storage/project_dasanthera_novaseq/results/trees/astral_CDS_noannotations.tre")

#if tip labels in species tree do not match those in gene trees...
for (i in 1:length(speciestree$tip.label)){
  speciestree$tip.label[i] <-
    str_split(str_split(speciestree$tip.label[i], "CDS_bws_")[[1]][2], ".fixed.fa")[[1]][1]
}

speciestree <- root(speciestree, outgroup = "mon_61-7_S440")
plot(speciestree)

#read in the set of gene trees of interest
intrees <- read.tree("~/project storage/project_dasanthera_novaseq/results/trees/combined_10kbwindowtrees.tre")


#read in the treepath names
treename_matrix <- read.table("~/project storage/project_dasanthera_novaseq/results/treemetrics/numbered_10kbtreepaths.txt")
colnames(treename_matrix) <- c('number','oldnames')


#define function to get information about the scaffold and region each window is from.
#(always give your files descriptive names :) )
addinfo <- function(j, treematrix, addtype) {
  if(addtype == "scaffold"){
    return(paste("scaffold_", 
                 strsplit((strsplit(treename_matrix[j,2], 'scaf_')[[1]][2]), '/')[[1]][1],
                 sep = '')
    )
  }
  if(addtype == "bpstart"){
    return(strsplit((strsplit(treematrix[j,2], 'bp_')[[1]][2]), '-')[[1]][1])
  }
  if(addtype == "bpend"){
    return(strsplit((strsplit(treename_matrix[j,2], '-')[[1]][2]), '.fa.treefile')[[1]][1])
  }
}


#add new columns to the treename matrix with this parsed information
treename_matrix$scaffolds <- sapply((1:nrow(treename_matrix)), addinfo,
                                    treematrix = treename_matrix,
                                    addtype = "scaffold")

treename_matrix$bpstart <- sapply((1:nrow(treename_matrix)), addinfo,
                                  treematrix = treename_matrix,
                                  addtype = "bpstart")

treename_matrix$bpend <- sapply((1:nrow(treename_matrix)), addinfo,
                                treematrix = treename_matrix,
                                addtype = "bpend")


#find RF distance for each tree. First, root trees to the outgroup.
for (i in 1:length(intrees)){
  intrees[[i]] <- root(phy = intrees[[i]], outgroup = "mon_61-7_S440")
}

#Then, estimate normalized RF distance and add to treename_matrix
treename_matrix$RFdist <- RF.dist(tree1 = intrees, tree2 = speciestree, normalize = T, rooted = T)

#write this to a file so we can load in again later
write.csv(treename_matrix, "RFdistance_10kbtrees_astralCDStree.csv", row.names = F)
##################################################



#Plotting RF distance in 10kb stretches along scaffold
#and overlay that onto plot of genic content
##################################################
#read in RF values and annofile info
RFfile <- read.csv("~/project storage/project_dasanthera_novaseq/results/treemetrics/RFdistance_10kbtrees_astral10kbtree.csv")
colnames(RFfile)[3] <- "chromosome"


annofile <- read.table("~/project storage/project_comparative_genome/genomes/davidsonii/annot_Pdavidsonii_1mb.gffread.genes.bed")
annofile <- annofile[,c(1,2,3,6,4)]
colnames(annofile) <- c("chromosome", "start", "end", "strand", "gene_id")

#make list of bad scaffolds
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")


#next, add midpoint info to the annotation file and treename file.
#this will be used to determine the number of CDS withiin tree windows.
annofile <- annofile %>%
  filter(!chromosome %in% badscafs) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(start, end)))/1000000) %>%
  ungroup()

RFfile <- RFfile %>%
  filter(!chromosome %in% badscafs) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(bpstart, bpend)))/1000000) %>%
  ungroup()


pdf("RFdist_genic_content_astral10kbtree.pdf", width = 12, height = 5)
ggplot() +
  geom_histogram(data = annofile, mapping = aes(x = midpoint, y = after_stat(ncount)),
                 binwidth = 1, alpha = 0.5) +
  geom_smooth(data = RFfile, mapping = aes(x = midpoint, y = RFdist),
              se = F, na.rm = T, linewidth = 0.5, method = "loess") +
  facet_grid(. ~ chromosome,
             scales = "free_x", space = "free_x") +
  xlab("Position on scaffold (Mb)") +
  scale_y_continuous(name = "10kb window tree normalized Robinson-Foulds distance",
                     sec.axis = sec_axis(~., name = "proportional gene count")) +
  ggtitle("loess-smoothed normalized RF distance (line) and binned gene content (hist)") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 8))
dev.off()

#loess, lm, glm, gam



