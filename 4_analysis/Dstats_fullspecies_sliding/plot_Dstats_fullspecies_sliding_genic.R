library(tidyverse)
#library(cowplot)
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_fullspecies_sliding_results/")


#want to generate a plot of f_d o y axis, gene density on x axis.
#read in the file, which was prepared with shell script "generate_dstat_genic_plotfiles.sh"

#get list of plotting files in the directory
filelist <- list.files(pattern = "plottingfile")

#set up df 
combdf <- data.frame()
for (i in 1:length(filelist)){
  tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
    mutate(test = str_split(filelist[i], "_localFstats")[[1]][1])
  
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
pdf("fdm_vs_genic_fraction.pdf")
ggplot(finaldata, aes(x = genic_fraction, y = absf_dM, group = speciescomp_fdM)) +
  geom_point(size = 0.2, alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE, col = "black", aes(group = comptype_fdM)) +
  geom_smooth(method = "lm", se = T, aes(col = speciescomp_fdM)) +
  facet_grid(~comptype_fdM) +
  ggtitle("f_dM vs. genic fraction") +
  scale_color_manual(values = colscheme)
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

write.csv(rbind(lm_output, eco_output), file = "lm_output.csv", row.names = F)

#value of 1 indicates a perfect fit
# value of 0 indicates model explains no variation in response variable



