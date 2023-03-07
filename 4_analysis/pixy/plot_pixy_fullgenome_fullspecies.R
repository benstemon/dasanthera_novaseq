setwd('~/Desktop/')
library(tidyverse)

#prepare list of regions containing genes of interest
genelist <- read.csv("~/project storage/project_dasanthera_novaseq/dasanthera_genelist.csv")

#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")




#PIXY 50KB
########################################################################
#set base directory for files
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_fullgenome_fullspecies_50kb/")

#add the genic fraction file and filter/mutate as needed
genicfractionfile <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_50kbwin_50kbslide.bed",
                                header = F, sep = " ",
                                col.names = c("chromosome", "window_pos_1", "window_pos_2", "genic_fraction"))
genicfractionfile <- genicfractionfile %>%
  mutate(window_pos_1 = window_pos_1+1)


#read in data -- pi, dxy, and fst
pidata <- read.delim("fullgenome_fullspecies_50kb_pi.txt", header = T)
dxydata <- read.delim("fullgenome_fullspecies_50kb_dxy.txt", header = T)
fstdata <- read.delim("fullgenome_fullspecies_50kb_fst.txt", header = T)


#join pidata with genic fraction, and fst+dxy with genic fraction
#also filter
plotting_pidata <- full_join(pidata, genicfractionfile) %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!pop %in% c("P_montanus", "P_lyallii")) %>%
  drop_na(c(avg_pi, genic_fraction)) %>%
  mutate(pop = as.character(str_sub_all(pop, 3, 5))) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  mutate(comptype = case_when(pop %in% c("car", "dav", "fru") ~ "bird",
                              pop %in% c("rup", "new") ~ "bee")) %>%
  ungroup()


#set up combinations for more easily categorizing interaction types
bird_bird = c("new_rup")
bee_bee = c("car_dav", "car_fru", "dav_fru")
bee_bird = c("dav_rup", "fru_rup", "dav_new", "fru_new", "car_rup", "car_new")

#do this for dxy-fst data
plotting_dxy_fstdata <- full_join(dxydata, fstdata) %>%
  full_join(., genicfractionfile) %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!pop1 %in% c("P_montanus", "P_lyallii")) %>%
  filter(!pop2 %in% c("P_montanus", "P_lyallii")) %>%
  drop_na(genic_fraction) %>%
  subset(avg_wc_fst >= 0) %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  na.omit() %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  mutate(inter = paste(str_sort(c(pop1, pop2))[1], str_sort(c(pop1, pop2))[2], sep = '_')) %>%
  mutate(comptype = case_when(inter %in% bird_bird ~ "bird_bird",
                              inter %in% bee_bee ~ "bee_bee",
                              inter %in% bee_bird ~ "bee_bird")) %>%
  ungroup()

#make a long version too
longplot_dxy_fst <- plotting_dxy_fstdata %>%
  pivot_longer(cols = c(avg_wc_fst, avg_dxy))


#generate an output matrix for variance, overall values, etc.
#for dxy, fst
write.csv(plotting_dxy_fstdata %>%
            group_by(inter) %>%
            summarize(fst_mean = mean(avg_wc_fst),
                      fst_min = min(avg_wc_fst),
                      fst_max = max(avg_wc_fst),
                      fst_sd = sd(avg_wc_fst),
                      dxy_mean = mean(avg_dxy),
                      dxy_min = min(avg_dxy),
                      dxy_max = max(avg_dxy),
                      dxy_sd = sd(avg_dxy)),
          "sumstats_dxy-fst_50kb.csv")

#write for pi
write.csv(plotting_pidata %>%
            group_by(pop) %>%
            summarize(pi_mean = mean(avg_pi),
                      pi_min = min(avg_pi),
                      pi_max = max(avg_pi),
                      pi_sd = sd(avg_pi)),
          "sumstats_pi_50kb.csv")


#generate boxplots for genome-wide metrics
a <- ggplot(longplot_dxy_fst, aes(x = inter, y = value)) +
  geom_boxplot(aes(col = comptype)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  scale_color_manual(values = c("bee_bee" = "blue",
                       "bee_bird" = "purple",
                       "bird_bird" = "red")) + 
  facet_grid(name ~.,
             scales = c("free_y"))

b <- ggplot(plotting_pidata, aes(x = pop, y = avg_pi)) +
  geom_boxplot(aes(col = comptype)) +
  scale_color_manual(values = c("bee" = "blue",
                                "bird" = "red")) +
  theme_bw()

library(cowplot)
pdf("BOXPLOT_genome-wide_50kb.pdf")
plot_grid(a, b,
          ncol = 1,
          rel_heights = c(1, 0.66))
dev.off()

#generate plots:
#1. pi vs. genic content vs genic content (per species average)
pdf("genic_fraction_vs_pi-loess_50kb.pdf", height =3.5)
ggplot(plotting_pidata, aes(x = genic_fraction, y = avg_pi, group = pop)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "loess", aes(col = pop)) +
  facet_grid(~pop)
dev.off()


#2. pi, fst, and dxy along the genome
pdf("pi_dxy_fst_along_genome_50kb.pdf", width = 10)
print(
  ggplot(plotting_pidata, aes(x = midpoint, y = avg_pi)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(pop ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)

print(
  ggplot(plotting_dxy_fstdata, aes(x = midpoint, y = avg_dxy)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(inter ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)

print(
  ggplot(plotting_dxy_fstdata, aes(x = midpoint, y = avg_wc_fst)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(inter ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)
dev.off()

#3. dxy and fst vs. genic content
#if png, use hw = 960, pointsize = 48
pdf("genic_fraction_vs_dxy-fst-loess_50kb.pdf")
ggplot(longplot_dxy_fst, aes(x = genic_fraction, y = value, group = inter)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "loess", aes(col = inter)) +
  facet_grid(name~comptype,
             scales = "free_y") +
  theme_bw()
dev.off()


#funsies. fst vs. dxy
pdf("dxy_vs_fst_50kb.pdf")
ggplot(plotting_dxy_fstdata, aes(x = avg_dxy, y = avg_wc_fst, group = inter)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "lm", aes(col = inter)) +
  facet_grid(~comptype)
dev.off()
########################################################################



#PIXY 10KB
########################################################################
#set base directory for files
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_fullgenome_fullspecies_10kb/")

#add the genic fraction file and filter/mutate as needed
genicfractionfile <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_10kbwin_10kbslide.bed",
                                header = F, sep = " ",
                                col.names = c("chromosome", "window_pos_1", "window_pos_2", "genic_fraction"))
genicfractionfile <- genicfractionfile %>%
  mutate(window_pos_1 = window_pos_1+1)


#read in data -- pi, dxy, and fst
pidata <- read.delim("fullgenome_fullspecies_10kb_pi.txt", header = T)
dxydata <- read.delim("fullgenome_fullspecies_10kb_dxy.txt", header = T)
fstdata <- read.delim("fullgenome_fullspecies_10kb_fst.txt", header = T)


#join pidata with genic fraction, and fst+dxy with genic fraction
#also filter
plotting_pidata <- full_join(pidata, genicfractionfile) %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!pop %in% c("P_montanus", "P_lyallii")) %>%
  drop_na(c(avg_pi, genic_fraction)) %>%
  mutate(pop = as.character(str_sub_all(pop, 3, 5))) %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  mutate(comptype = case_when(pop %in% c("car", "dav", "fru") ~ "bird",
                              pop %in% c("rup", "new") ~ "bee")) %>%
  ungroup()


#set up combinations for more easily categorizing interaction types
bird_bird = c("new_rup")
bee_bee = c("car_dav", "car_fru", "dav_fru")
bee_bird = c("dav_rup", "fru_rup", "dav_new", "fru_new", "car_rup", "car_new")

#do this for dxy-fst data
plotting_dxy_fstdata <- full_join(dxydata, fstdata) %>%
  full_join(., genicfractionfile) %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!pop1 %in% c("P_montanus", "P_lyallii")) %>%
  filter(!pop2 %in% c("P_montanus", "P_lyallii")) %>%
  drop_na(genic_fraction) %>%
  subset(avg_wc_fst >= 0) %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  na.omit() %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(window_pos_1, window_pos_2)))/1000000) %>%
  mutate(inter = paste(str_sort(c(pop1, pop2))[1], str_sort(c(pop1, pop2))[2], sep = '_')) %>%
  mutate(comptype = case_when(inter %in% bird_bird ~ "bird_bird",
                              inter %in% bee_bee ~ "bee_bee",
                              inter %in% bee_bird ~ "bee_bird")) %>%
  ungroup()

#make a long version too
longplot_dxy_fst <- plotting_dxy_fstdata %>%
  pivot_longer(cols = c(avg_wc_fst, avg_dxy))


#generate an output matrix for variance, overall values, etc.
#for dxy, fst
write.csv(plotting_dxy_fstdata %>%
            group_by(inter) %>%
            summarize(fst_mean = mean(avg_wc_fst),
                      fst_min = min(avg_wc_fst),
                      fst_max = max(avg_wc_fst),
                      fst_sd = sd(avg_wc_fst),
                      dxy_mean = mean(avg_dxy),
                      dxy_min = min(avg_dxy),
                      dxy_max = max(avg_dxy),
                      dxy_sd = sd(avg_dxy)),
          "sumstats_dxy-fst_10kb.csv")

#write for pi
write.csv(plotting_pidata %>%
            group_by(pop) %>%
            summarize(pi_mean = mean(avg_pi),
                      pi_min = min(avg_pi),
                      pi_max = max(avg_pi),
                      pi_sd = sd(avg_pi)),
          "sumstats_pi_10kb.csv")


#generate boxplots for genome-wide metrics
a <- ggplot(longplot_dxy_fst, aes(x = inter, y = value)) +
  geom_boxplot(aes(col = comptype)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  scale_color_manual(values = c("bee_bee" = "blue",
                                "bee_bird" = "purple",
                                "bird_bird" = "red")) + 
  facet_grid(name ~.,
             scales = c("free_y"))

b <- ggplot(plotting_pidata, aes(x = pop, y = avg_pi)) +
  geom_boxplot(aes(col = comptype)) +
  scale_color_manual(values = c("bee" = "blue",
                                "bird" = "red")) +
  theme_bw()

library(cowplot)
pdf("BOXPLOT_genome-wide_10kb.pdf")
plot_grid(a, b,
          ncol = 1,
          rel_heights = c(1, 0.66))
dev.off()

#generate plots:
#1. pi vs. genic content vs genic content (per species average)
pdf("genic_fraction_vs_pi-loess_10kb.pdf", height =3.5)
ggplot(plotting_pidata, aes(x = genic_fraction, y = avg_pi, group = pop)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "loess", aes(col = pop)) +
  facet_grid(~pop)
dev.off()


#2. pi, fst, and dxy along the genome
pdf("pi_dxy_fst_along_genome_10kb.pdf", width = 10)
print(
  ggplot(plotting_pidata, aes(x = midpoint, y = avg_pi)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(pop ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)

print(
  ggplot(plotting_dxy_fstdata, aes(x = midpoint, y = avg_dxy)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(inter ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)

print(
  ggplot(plotting_dxy_fstdata, aes(x = midpoint, y = avg_wc_fst)) +
    geom_point(size = 0.2, alpha = 0.4) +
    geom_smooth(method = "loess", se = F) +
    facet_grid(inter ~ chromosome,
               scales = "free_x", space = "free_x") +
    xlab("Position on scaffold (Mb)") +
    theme_bw() +
    theme(strip.text.x = element_text(size = 7.5))
)
dev.off()

#3. dxy and fst vs. genic content
#if png, use hw = 960, pointsize = 48
pdf("genic_fraction_vs_dxy-fst-loess_10kb.pdf")
ggplot(longplot_dxy_fst, aes(x = genic_fraction, y = value, group = inter)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "loess", aes(col = inter)) +
  facet_grid(name~comptype,
             scales = "free_y") +
  theme_bw()
dev.off()


#funsies. fst vs. dxy
pdf("dxy_vs_fst_50kb.pdf")
ggplot(plotting_dxy_fstdata, aes(x = avg_dxy, y = avg_wc_fst, group = inter)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "lm", aes(col = inter)) +
  facet_grid(~comptype)
dev.off()
########################################################################




#keep this, as it shows how to incorporate genes of interest
############################################################
#read in raw data
rawdata <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixyout_fullgenome_fullspecies_50kb/fullgenome_fullspecies_50kb_dxy.txt",
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
png("dxy_intermorph_50kb.png", width = 3000, height = 3000, res = 300)
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




#plot dxy -- intramorph
png("dxy_intramorph_50kb.png", width = 3000, height = 3000, res = 300)
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


#plot genome-wide values for dxy in density plots
t <- rbind(sub_intermorph, sub_intramorph)
pdf("genome-wide_dxy_50kb.pdf")
ggplot(t, aes(x = avg_dxy)) +
  geom_density(aes(group = inter, col = inter))
dev.off()


############################################################
gtable_show_layout(gt)
gt$heights[8]

