setwd('~/Desktop/')
library(tidyverse)


#prepare list of regions containing genes of interest
#genelist <- read.csv("~/project storage/project_dasanthera_novaseq/dasanthera_genelist.csv")


#read in 10kb data
rawdata <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixyout_fullgenome_fullspecies_10kb/fullgenome_fullspecies_10kb_dxy.txt",
                      header = T)

#do some basic cleanup and remove missing data
#first generate objects to replace poorly chosen species names
#and bad scaffolds
badtaxa <- c("P_montanus", "P_lyallii")
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")



#names sorted alphabetically already so no need to sort for speciescomp
#also filtering only for tests that include fruticosus
fdata <- rawdata %>%
  na.omit() %>%
  filter(!pop1 %in% badtaxa) %>%
  filter(!pop2 %in% badtaxa) %>%
  filter(!chromosome %in% badscafs) %>%
  mutate_at(c("pop1", "pop2"), ~str_replace(., "P_cardwellii", "car")) %>%
  mutate_at(c("pop1", "pop2"), ~str_replace(., "P_rupicola", "rup")) %>%
  mutate_at(c("pop1", "pop2"), ~str_replace(., "P_newberryi", "new")) %>%
  mutate_at(c("pop1", "pop2"), ~str_replace(., "P_davidsonii", "dav")) %>%
  mutate_at(c("pop1", "pop2"), ~str_replace(., "P_fruticosus", "fru")) %>%
  mutate(speciescomp = paste(pop1, pop2, sep = "_")) %>%
  filter(., grepl("fru", speciescomp)) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  ungroup()



#plot differences 
pdf("dxy_frut_multicomp1.pdf")
ggplot(fdata, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.02, alpha = 0.4) +
  geom_smooth(method = "loess") +
  facet_grid(speciescomp ~ chromosome,
             scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = "5")) +
  xlab("Position on chromosome (Mb)")
dev.off()

pdf("dxy_frut_multicomp2.pdf", height = 3)
ggplot(fdata, aes(x = midpoint, y = avg_dxy, col = speciescomp)) +
  geom_smooth(method = "loess", linewidth = 0.5) +
  facet_grid(~chromosome,
             scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size = "5")) +
  xlab("Position on chromosome (Mb)")
dev.off()

#10kb PIXY dxy
############################################################
#read in raw data
rawdata <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixyout_fullgenome_fullspecies_10kb/fullgenome_fullspecies_10kb_dxy.txt",
                      header = T)

#make additional colum to specify the species interaction
test <- rawdata %>%
  na.omit() %>%
  rowwise() %>%
  mutate(inter = paste(str_sort(c(pop1, pop2))[1], str_sort(c(pop1, pop2))[2], sep = '_x_')) %>%
  ungroup()


#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")

#set up subset data with only cross-morph comparisons
intermorph = c("P_fruticosus_x_P_rupicola",
               "P_davidsonii_x_P_rupicola",
               "P_davidsonii_x_P_newberryi",
               "P_cardwellii_x_P_rupicola",
               "P_cardwellii_x_P_newberryi",
               "P_fruticosus_x_P_newberryi")
#also filter for low snps and bad chromosomes, and bad fst estimates
sub_intermorph <- test %>%
  filter(no_sites > 100) %>%
  filter(inter %in% intermorph) %>%
  filter(!chromosome %in% badscafs) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  ungroup()

#make a sub_df which identifies which regions in the cross-morph comps contain genes of interest
genepoints_intermorph <- sub_intermorph %>%
  rowwise() %>%
  mutate(genename = case_when((between(genelist$midpoint[1], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[1]) ~ genelist$gene[1],
                              (between(genelist$midpoint[2], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[2]) ~ genelist$gene[2],
                              (between(genelist$midpoint[3], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[3]) ~ genelist$gene[3],
                              (between(genelist$midpoint[4], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[4]) ~ genelist$gene[4],
                              (between(genelist$midpoint[5], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[5]) ~ genelist$gene[5])) %>%
  ungroup() %>%
  na.omit()



#set up subset data with only within-morph comparisons
intramorph = c("P_cardwellii_x_P_davidsonii",
               "P_cardwellii_x_P_fruticosus",
               "P_newberryi_x_P_rupicola",
               "P_davidsonii_x_P_fruticosus")
sub_intramorph <- test %>%
  filter(no_sites > 100) %>%
  filter(inter %in% intramorph) %>%
  filter(!chromosome %in% badscafs) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  ungroup()

#make the genes of interest sub_df for this as well
genepoints_intramorph <- sub_intramorph %>%
  rowwise() %>%
  mutate(genename = case_when((between(genelist$midpoint[1], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[1]) ~ genelist$gene[1],
                              (between(genelist$midpoint[2], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[2]) ~ genelist$gene[2],
                              (between(genelist$midpoint[3], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[3]) ~ genelist$gene[3],
                              (between(genelist$midpoint[4], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[4]) ~ genelist$gene[4],
                              (between(genelist$midpoint[5], window_pos_1, window_pos_2) & 
                                 chromosome == genelist$scaffold[5]) ~ genelist$gene[5])) %>%
  ungroup() %>%
  na.omit()


####PLOTTING
library(grid)
library(gtable)

#plot dxy -- inter-morph
png("dxy_intermorph_10kb.png", width = 3000, height = 3000, res = 300)
a <- ggplot(sub_intermorph, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(inter ~ chromosome,
             scales = "free_x", space = "free_x") +
  theme(strip.text.y = element_text(angle = 0)) +
  xlab("Position on scaffold (Mb)") +
  geom_point(data = genepoints_intermorph, aes(x = midpoint, y = avg_dxy, colour = genename))

#alter some of the plot values so they compare better
gt = ggplot_gtable(ggplot_build(a))
gt$widths[20] = grid::unit(3, "cm")
gt$heights[c(8,10,12,14,16,18)] = grid::unit(3.5, "cm")
grid.draw(gt)
dev.off()


#plot genome-wide values for dxy in density plots
t <- rbind(sub_intermorph, sub_intramorph)
pdf("genome-wide_dxy_10kb.pdf")
ggplot(t, aes(x = avg_dxy)) +
  geom_density(aes(group = inter, col = inter))
dev.off()



#plot dxy -- intramorph
png("dxy_intramorph_10kb.png", width = 3000, height = 3000, res = 300)
b <- ggplot(sub_intramorph, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(inter ~ chromosome,
             scales = "free_x", space = "free_x") +
  theme(strip.text.y = element_text(angle = 0)) +
  xlab("Position on scaffold (Mb)") +
  geom_point(data = genepoints_intramorph, aes(x = midpoint, y = avg_dxy, colour = genename))

#alter some of the plot values so they compare better
gt = ggplot_gtable(ggplot_build(b))
gt$widths[20] = grid::unit(3, "cm")
gt$heights[c(8,10,12,14)] = grid::unit(3.5, "cm")
grid.draw(gt)

dev.off()


############################################################
gtable_show_layout(gt)
gt$heights[8]