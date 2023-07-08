setwd('~/Desktop/')
library(tidyverse)

#set base directory for files
setwd("~/project storage/project_dasanthera_novaseq/results/dxy_outliers_hybridzone_results/")

#prepare list of regions containing genes of interest
genelist <- read.csv("~/project storage/project_dasanthera_novaseq/dasanthera_genelist.csv")

#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")

#add the genic fraction file and filter/mutate as needed
genicfractionfile <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_10kbwin_10kbslide.bed",
                                header = F, sep = " ",
                                col.names = c("chromosome", "window_pos_1", "window_pos_2", "genic_fraction")) %>%
  mutate(window_pos_1 = window_pos_1+1)
  


#read in data of interest
dxydata <- read.delim("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_fullgenome_fullspecies_10kb/fullgenome_fullspecies_10kb_dxy.txt", header = T)




##data processing
#filter data to include pops of interest, bad scafs, etc.
maindata <- left_join(dxydata, genicfractionfile,
                      by = c("chromosome", "window_pos_1", "window_pos_2")) %>%
  filter(!chromosome %in% badscafs) %>%
  na.omit() %>%
  mutate(midpoint = ceiling((window_pos_1+window_pos_2)/2)/1000000,
         pop1 = substr(pop1, 3, 5),
         pop2 = substr(pop2, 3, 5),
         chromosome = gsub("fold", "", chromosome),
         comp = paste(pop1, pop2, sep="_"),
         comp = ifelse(comp %in% c("new_dav", "dav_new"), "dav_new", comp),
         comp = ifelse(comp %in% c("new_rup", "rup_new"), "new_rup", comp),
         comp = ifelse(comp %in% c("rup_dav", "dav_rup"), "dav_rup", comp),
         comp = ifelse(comp %in% c("rup_fru", "fru_rup"), "fru_rup", comp)) %>%
  filter(comp %in% c("dav_new", "new_rup", "dav_rup", "fru_rup")) %>%
  group_by(comp) %>%
  mutate(zscore_unfiltered = (avg_dxy - mean(avg_dxy)) / sd(avg_dxy))



####
#find mean dxy values for each comparison
mean_dxy_unfiltered <- maindata %>%
  group_by(comp) %>%
  summarize(mean_avg_dxy = mean(avg_dxy, na.rm = TRUE))
write.csv(mean_dxy_unfiltered, file = "mean_dxy_unfiltered_hybridzonepops.csv")

#find subset of outlier loci
outliers_unfiltered <- maindata %>% filter(zscore_unfiltered > 4)
write.csv(outliers_unfiltered, "4z_dxy_hybridzone_outliers_unfiltered.csv")

#output format usable for .bed
write_delim(outliers_unfiltered[,c(3:5, 13)] %>% mutate(chromosome = gsub("scaf", "scaffold", chromosome)),
            "4z_dxy_hybridzone_outliers_unfiltered.bed",
            delim = '\t', col_names = F)


#plot raw manhattan plots for each comparison
png("dxy_hybridzonepops_4z_outliers_unfiltered.png", width = 10, height = 4.5, units = "in", res = 400)
ggplot(maindata, aes(x = midpoint, y = avg_dxy)) +
  facet_grid(comp~chromosome, space = "free_x", scales = "free_x") +
  geom_point(size = 0.3) +
  geom_point(data = outliers_unfiltered, color = "red", size = 0.35) +
  theme_bw() +
  theme(panel.spacing.x = unit(0.01, "in")) +
  xlab("Position on scaffold (Mb)") +
  ylab(expression(paste("Average ", italic(d[xy]))))
dev.off()

#saveRDS(maindata, file="~/project storage/project_dasanthera_novaseq/plot-making-compendium/unfiltered_dxy_hybridzonepops.obj")
#saveRDS(outliers_unfiltered, file="~/project storage/project_dasanthera_novaseq/plot-making-compendium/unfiltered_dxy_hybridzonepops_outliers.obj")
####




#filter data more thoroughly for missingness, etc.
#1. only sites with fewer missing counts than count comparisons
filterdata <- maindata %>%
  filter(count_missing < count_comparisons) %>%
  group_by(comp) %>%
  mutate(zscore_filtered = (avg_dxy - mean(avg_dxy)) / sd(avg_dxy))

####
#find mean dxy values for each comparison
mean_dxy_filtered <- filterdata %>%
  group_by(comp) %>%
  summarize(mean_avg_dxy = mean(avg_dxy, na.rm = TRUE))
write.csv(mean_dxy_filtered, file = "mean_dxy_filtered_hybridzonepops.csv")

#find subset of outlier loci
outliers_filtered <- filterdata %>% filter(zscore_filtered > 4)
write.csv(outliers_filtered, "4z_dxy_hybridzone_outliers_filtered.csv")

#filter for duplicates and output format usable for .bed
for (i in unique(outliers_filtered$comp)){
  write_delim(outliers_filtered[,c(3:5,13)] %>%
                mutate(chromosome = gsub("scaf", "scaffold", chromosome)) %>%
                filter(comp %in% i),
              paste("4z_dxy_hybridzone_outliers_filtered_", i, ".bed", sep=""),
              delim = "\t", col_names = F)
}

#plot raw manhattan plots for each comparison
png("dxy_hybridzonepops_4z_outliers_filtered.png", width = 10, height = 4.5, units = "in", res = 400)
ggplot(filterdata, aes(x = midpoint, y = avg_dxy)) +
  facet_grid(comp~chromosome, space = "free_x", scales = "free_x") +
  geom_point(size = 0.3) +
  geom_point(data = outliers_filtered, color = "red", size = 0.35) +
  theme_bw() +
  theme(panel.spacing.x = unit(0.01, "in")) +
  xlab("Position on scaffold (Mb)") +
  ylab(expression(paste("Average ", italic(d[xy]))))
dev.off()


#regression for dxy vs. genic fraction
#fit linear models
models <- filterdata %>%
  group_by(comp) %>%
  do(model = lm(avg_dxy ~ genic_fraction, data = .))

#make table with relevant information about linear models
lm_output <- data.frame(matrix(nrow = length(models$comp), ncol = 4))
colnames(lm_output) <- c("comp", "coeff", "r2", "pval")

for (i in 1: length(models$comp)){
  lm_output[i,1:4] <- c(models$comp[i],
                        summary(models$model[[i]])$coefficients[2],
                        summary(models$model[[i]])$r.squared,
                        anova(models$model[[i]])$`Pr(>F)`[1])
}
write.csv(lm_output, file = "regression_filtered_hybridzonepops_linearmodels_dxy_vs_genicfraction.csv")


#generate plot for this
png("dxy_hybridzonepops_vs._genicfraction_filtered.png", width = 2.2, height = 5.3, units = "in", res = 400)
ggplot(filterdata, aes(x = genic_fraction, y = avg_dxy)) +
  geom_point(size = 0.3, color = "gray40", alpha = 0.5) +
  geom_smooth(method = "lm", se = F, color = "firebrick") +
  facet_grid(comp~.) +
  theme_bw() +
  theme(panel.spacing.x = unit(0.01, "in")) +
  xlab("Genic fraction") +
  ylab(expression(paste("Average ", italic(d[xy]))))
dev.off()

saveRDS(filterdata, file="~/project storage/project_dasanthera_novaseq/plot-making-compendium/filtered_dxy_hybridzonepops.obj")
saveRDS(outliers_filtered, file="~/project storage/project_dasanthera_novaseq/plot-making-compendium/filtered_dxy_hybridzonepops_outliers.obj")
####

