library(tidyverse)
library(ggrepel)

setwd("~/project storage/project_dasanthera_novaseq/results/PCA_genome/")


# read in data
pca <- read_table("PCA.eigenvec", col_names = FALSE)
eigenval <- scan("PCA.eigenval")


#clean up data
pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#add species and pop columns, with some individual-specific criteria
pca <- pca %>%
  rowwise() %>%
  mutate(pop = str_split(str_split(ind, 'bws_')[[1]][2], '-')[[1]][1]) %>%
  mutate(species = case_when(pop == "dav_118" ~ "fru",
                             pop != "dav_118" ~ str_split(pop, '_')[[1]][1]))


# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)



# plot pca
a <- ggplot(pca, aes(x = PC1, y = PC2, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("genome-wide PCA -- plink version -- PC1 and 2")

#PC2 and 3
b <- ggplot(pca, aes(x = PC2, y = PC3, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ggtitle("genome-wide PCA -- plink version -- PC2 and 3")

#PC3 and PC4
c <- ggplot(pca, aes(x = PC3, y = PC4, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) +
  ggtitle("genome-wide PCA -- plink version -- PC3 and 4")

pdf("genome_PCA-plink.pdf")
print(a)
print(b)
print(c)
dev.off()



##### Try loading in the vcf and making PCA in R --
#the plink one is a little weird
setwd("~/project storage/project_dasanthera_novaseq/results/PCA_genome/")
library(tidyverse)
library(adegenet)
library(vcfR)
library(ggrepel)
library(gridExtra)

#read in VCF and convert to adegenet genlight object
invcf <- read.vcfR("PCA-ready_biallelic_ld-filtered.vcf")
my_genlight <- vcfR2genlight(invcf)


#make df to then use for adding ind and species information
#keep in mind special case for dav_118
t <- data.frame(ind = my_genlight$ind.names) %>%
  rowwise() %>%
  mutate(pop = str_split(str_split(ind, 'bws_')[[1]][2], '-')[[1]][1]) %>%
  mutate(species = case_when(pop == "dav_118" ~ "fru",
                             pop != "dav_118" ~ str_split(pop, '_')[[1]][1]))

#replace or update pertinent genlight attributes
my_genlight$ind.names <- as.vector(t$pop)


#make pca
pca <- glPca(my_genlight, nf = 18)

#check that eigenvalues are not negative
pca$eig
sum(pca$eig)
pve <- data.frame(PC = 1:length(pca$eig), pve = pca$eig/sum(pca$eig)*100)


#need to transform this a bit because the the PCA object isn't in the best format
#so just add PCA values to the t object
r <- cbind(t, pca$scores[,1], pca$scores[,2], pca$scores[,3], pca$scores[,4])
colnames(r)[4:7] <- c('PC1', 'PC2', 'PC3', 'PC4')


#plot again
#PC1 and 2
a <- ggplot(r, aes(x = PC1, y = PC2, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ggtitle("genome-wide PCA -- adegenet version -- PC1 and 2")

#PC2 and 3
b <- ggplot(r, aes(x = PC2, y = PC3, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ggtitle("genome-wide PCA -- adegenet version -- PC2 and 3")

#PC3 and 4
c <- ggplot(r, aes(x = PC3, y = PC4, col = species)) + geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) +
  ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)")) +
  ggtitle("genome-wide PCA -- adegenet version -- PC3 and 4")

pdf("genome_PCA-adegenet.pdf")
print(a)
print(b)
print(c)
dev.off()


