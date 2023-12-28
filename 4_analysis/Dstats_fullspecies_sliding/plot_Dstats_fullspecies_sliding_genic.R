library(tidyverse)
#library(cowplot)
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_fullspecies_sliding_results/v1")



#want to generate a plot of f_d on y axis, gene density on x axis.
#read in the file, which was prepared with shell script "generate_dstat_genic_plotfiles.sh"

#Version 1, with fru106 included:
############################################################

#get list of plotting files in the directory
filelist <- list.files(recursive = T, pattern = "10kb_10000_2500_plottingfile")

#set up df 
combdf <- data.frame()
for (i in 1:length(filelist)){
  tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
    mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
  
  combdf <- rbind(combdf, tmpdata)
}


# To just analyze each cross separately:
# this creates new statistic absf_dM, which is always positive.
# This is because fdM can be positive or negative.
# When positive, comparison is P2-P3
# When negative, comparison is P1-P3
finaldata <- combdf %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(windowStart, windowEnd)))) %>%
  mutate(absf_dM = abs(f_dM)) %>%
  mutate(absf_d = abs(f_d)) %>%
  mutate(speciescomp_fdM = case_when((f_dM < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                     (f_dM > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  mutate(speciescomp_fd = case_when((f_d < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                    (f_d > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  ungroup()


# Save this df to make a better plot later:
save(finaldata, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fdm_vs_genic.obj")

#plot results
pdf("allcomps_fdm_vs_genic.pdf", height = 6, width = 15)
ggplot(finaldata, aes(x = genic_fraction, y = absf_dM, group = speciescomp_fdM)) +
  geom_point(size = 0.6, alpha = 0.2) +
  geom_smooth(method = "lm", se = T, linewidth = 2) +
  facet_wrap(~speciescomp_fdM,nrow=2,scales="free") +
  ggtitle("f_dM vs. genic fraction") +
  ylim(0,0.25) +
  theme_bw() +
  theme(text = element_text(size = 24))
dev.off()



# To make categories and do it that way:
######################################################################
#make a category to bin comparisons based on introgression "type"
#be sure to list species comps in alphabetical order (this is how they will be sorted)
bird_bird = c("new_rup")
bee_bee = c("car_dav", "car_fru")
bee_bird = c("dav_rup", "fru_rup", "dav_new", "fru_new", "car_rup")

#need to add a bit more information, such as midpoint for plotting, 
#and a category for the statistic (if positive, P2-P3, if negative, P1-P3)
finaldata <- combdf %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(windowStart, windowEnd)))) %>%
  mutate(absf_dM = abs(f_dM)) %>%
  mutate(absf_d = abs(f_d)) %>%
  mutate(speciescomp_fdM = case_when((f_dM < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                 (f_dM > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  mutate(speciescomp_fd = case_when((f_d < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                     (f_d > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  mutate(comptype_fdM = case_when(speciescomp_fdM %in% bird_bird ~ "bird_bird",
                                  speciescomp_fdM %in% bee_bee ~ "bee_bee",
                                  speciescomp_fdM %in% bee_bird ~ "bee_bird")) %>%
  mutate(comptype_fd = case_when(speciescomp_fd %in% bird_bird ~ "bird_bird",
                                 speciescomp_fd %in% bee_bee ~ "bee_bee",
                                 speciescomp_fd %in% bee_bird ~ "bee_bird")) %>%
  ungroup()


#make color scheme for plotting
colscheme <- c("car_dav" = "#0000ff", "car_fru" = "#006400", "car_rup" = "#ff0000",
               "dav_new" = "#ffd700", "dav_rup" = "#00ff00", "fru_new" = "#00ffff",
               "fru_rup" = "#ff00ff", "new_rup" = "#ffb6c1")


#plot results
pdf("fdm_vs_genic_fraction_v1.pdf", height = 6, width = 15)
ggplot(finaldata, aes(x = genic_fraction, y = absf_dM, group = speciescomp_fdM)) +
  geom_point(size = 0.6, alpha = 0.2) +
  geom_smooth(method = "lm", se = T, aes(col = speciescomp_fdM), linewidth = 2) +
  facet_grid(~comptype_fdM) +
  ggtitle("f_dM vs. genic fraction") +
  scale_color_manual(values = colscheme) +
  ylim(0,0.25) +
  theme_bw() +
  theme(text = element_text(size = 24))
dev.off()
######################################################################

#fit linear models
models <- finaldata %>%
  group_by(speciescomp_fdM) %>%
  do(model = lm(absf_dM ~ genic_fraction, data = .))


#make table with relevant information about linear models
lm_output <- data.frame(matrix(nrow = length(models$speciescomp_fdM), ncol = 4))
colnames(lm_output) <- c("species", "coeff", "r2", "pval")

for (i in 1: length(models$speciescomp_fdM)){
  lm_output[i,1:4] <- c(models$speciescomp_fdM[i],
                        summary(models$model[[i]])$coefficients[2],
                        summary(models$model[[i]])$r.squared,
                        anova(models$model[[i]])$`Pr(>F)`[1])
}

#also want to add in the regressions at the type of interaction
# only if you split up by the eco type
ecomodels <- finaldata %>%
  group_by(comptype_fdM) %>%
  do(model = lm(absf_dM ~ genic_fraction, data = .))

#make df for this and add to it 
eco_output <- data.frame(matrix(nrow = length(ecomodels$comptype_fdM), ncol = 4))
colnames(eco_output) <- c("species", "coeff", "r2", "pval")

for (i in 1: length(ecomodels$comptype_fdM)){
  eco_output[i,1:4] <- c(ecomodels$comptype_fdM[i],
                        summary(ecomodels$model[[i]])$coefficients[2],
                        summary(ecomodels$model[[i]])$r.squared,
                        anova(ecomodels$model[[i]])$`Pr(>F)`[1])
}


#combine the two dfs together, and write to output

write.csv(rbind(lm_output, eco_output), file = "lm_output_v1.csv", row.names = F)
############################################################


#Version 2, with no fru106
############################################################

#get list of plotting files in the directory
filelist <- list.files(recursive = T, pattern = "nofru106_10000_2500_plottingfile")

#set up df 
combdf <- data.frame()
for (i in 1:length(filelist)){
  tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
    mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
  
  combdf <- rbind(combdf, tmpdata)
}


#make a category to bin comparisons based on introgression "type"
#be sure to list species comps in alphabetical order (this is how they will be sorted)
bird_bird = c("new_rup")
bee_bee = c("car_dav", "car_fru")
bee_bird = c("dav_rup", "fru_rup", "dav_new", "fru_new", "car_rup")

#need to add a bit more information, such as midpoint for plotting, 
#and a category for the statistic (if positive, P2-P3, if negative, P1-P3)
finaldata <- combdf %>%
  rowwise() %>%
  mutate(midpoint = ceiling(mean(c(windowStart, windowEnd)))) %>%
  mutate(absf_dM = abs(f_dM)) %>%
  mutate(absf_d = abs(f_d)) %>%
  mutate(speciescomp_fdM = case_when((f_dM < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                     (f_dM > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  mutate(speciescomp_fd = case_when((f_d < 0) ~ paste(sort(str_split_1(test, "_")[c(1,3)]), collapse = "_"),
                                    (f_d > 0) ~ paste(sort(str_split_1(test, "_")[c(2,3)]), collapse = "_"))) %>%
  mutate(comptype_fdM = case_when(speciescomp_fdM %in% bird_bird ~ "bird_bird",
                                  speciescomp_fdM %in% bee_bee ~ "bee_bee",
                                  speciescomp_fdM %in% bee_bird ~ "bee_bird")) %>%
  mutate(comptype_fd = case_when(speciescomp_fd %in% bird_bird ~ "bird_bird",
                                 speciescomp_fd %in% bee_bee ~ "bee_bee",
                                 speciescomp_fd %in% bee_bird ~ "bee_bird")) %>%
  ungroup()


#make color scheme for plotting
colscheme <- c("car_dav" = "#0000ff", "car_fru" = "#006400", "car_rup" = "#ff0000",
               "dav_new" = "#ffd700", "dav_rup" = "#00ff00", "fru_new" = "#00ffff",
               "fru_rup" = "#ff00ff", "new_rup" = "#ffb6c1")


#plot results
pdf("fdm_vs_genic_fraction_nofru106.pdf")
ggplot(finaldata, aes(x = genic_fraction, y = absf_dM, group = speciescomp_fdM)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE, col = "black", aes(group = comptype_fdM)) +
  geom_smooth(method = "lm", se = T, aes(col = speciescomp_fdM)) +
  facet_grid(~comptype_fdM) +
  ggtitle("f_dM vs. genic fraction") +
  scale_color_manual(values = colscheme) +
  ylim(0,0.3)
dev.off()


#fit linear models
models <- finaldata %>%
  group_by(speciescomp_fdM) %>%
  do(model = lm(absf_dM ~ genic_fraction, data = .))


#make table with relevant information about linear models
lm_output <- data.frame(matrix(nrow = length(models$speciescomp_fdM), ncol = 4))
colnames(lm_output) <- c("species", "coeff", "r2", "pval")

for (i in 1: length(models$speciescomp_fdM)){
  lm_output[i,1:4] <- c(models$speciescomp_fdM[i],
                        summary(models$model[[i]])$coefficients[2],
                        summary(models$model[[i]])$r.squared,
                        anova(models$model[[i]])$`Pr(>F)`[1])
}

#also want to add in the regressions at the type of interaction
ecomodels <- finaldata %>%
  group_by(comptype_fdM) %>%
  do(model = lm(absf_dM ~ genic_fraction, data = .))

#make df for this and add to it 
eco_output <- data.frame(matrix(nrow = length(ecomodels$comptype_fdM), ncol = 4))
colnames(eco_output) <- c("species", "coeff", "r2", "pval")

for (i in 1: length(ecomodels$comptype_fdM)){
  eco_output[i,1:4] <- c(ecomodels$comptype_fdM[i],
                         summary(ecomodels$model[[i]])$coefficients[2],
                         summary(ecomodels$model[[i]])$r.squared,
                         anova(ecomodels$model[[i]])$`Pr(>F)`[1])
}


#combine the two dfs together, and write to output

write.csv(rbind(lm_output, eco_output), file = "lm_output_nofru106.csv", row.names = F)
############################################################

