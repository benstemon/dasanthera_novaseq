setwd('~/Desktop/')
library(tidyverse)

#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")


#read in data
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_hybridzone_divergence/")

dxy_10kb <- read_delim("fullgenome_fullspecies_10kb_dxy.txt")
dxy_50kb <- read_delim("fullgenome_fullspecies_50kb_dxy.txt")

#filter data to make plots
filtered_10kb <- dxy_10kb %>%
  filter(!chromosome %in% badscafs) %>%
  na.omit() %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  mutate(midpoint = ceiling((window_pos_1-1 + window_pos_2)/2)/1000000) %>%
  filter(pop1 %in% c("dav", "new", "rup") & pop2 %in% c("dav", "new", "rup")) %>%
  mutate(pop_comparison = paste(pop1, pop2, sep = "_")) %>%
  filter(pop_comparison != "new_rup") %>%
  mutate(chromosome = gsub("fold", "", chromosome))

filtered_50kb <- dxy_50kb %>%
  filter(!chromosome %in% badscafs) %>%
  na.omit() %>%
  mutate(pop1 = as.character(str_sub_all(pop1, 3, 5))) %>%
  mutate(pop2 = as.character(str_sub_all(pop2, 3, 5))) %>%
  mutate(midpoint = ceiling((window_pos_1-1 + window_pos_2)/2)/1000000) %>%
  filter(pop1 %in% c("dav", "new", "rup") & pop2 %in% c("dav", "new", "rup")) %>%
  mutate(pop_comparison = paste(pop1, pop2, sep = "_")) %>%
  filter(pop_comparison != "new_rup") %>%
  mutate(chromosome = gsub("fold", "", chromosome))


#identify outlier dxy values based on z score
outliers_10kb <- filtered_10kb %>%
  mutate(zscore = ((avg_dxy - mean(avg_dxy))/sd(avg_dxy))) %>%
  filter(abs(zscore) > 3)

outliers_50kb <- filtered_50kb %>%
  mutate(zscore = ((avg_dxy - mean(avg_dxy))/sd(avg_dxy))) %>%
  filter(abs(zscore) > 3)


#plot these results
ggplot(filtered_10kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = outliers_10kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x')

ggplot(filtered_50kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = outliers_50kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x')



#check if any outlier windows are shared between the two comparisons
shared_outliers_10kb <- outliers_10kb %>%
  group_by(chromosome, midpoint) %>%
  filter(n() > 1) %>%
  ungroup()

shared_outliers_50kb <- outliers_50kb %>%
  group_by(chromosome, midpoint) %>%
  filter(n() > 1) %>%
  ungroup()


#plot shared outliers
a <- ggplot(filtered_10kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers_10kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  labs(x = "Position on chromosome (Mbp)",
       y = expression(paste("Average ", italic(d[xy])))) +
  theme(panel.spacing.x = unit(0.01, "in")) +
  ggtitle("Shared dxy outliers 10kb windows")
write.csv(shared_outliers_10kb, "shared_dxy_outliers_10kb.csv", row.names = F)
png(filename = "shared_dxy_outliers_10kb.png",
    height = 2.5, width = 9, units = "in", res = 400)
a
dev.off()

b <- ggplot(filtered_50kb, aes(x = midpoint, y = avg_dxy)) +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers_50kb, colour = "red", size = 0.4) +
  facet_grid(pop_comparison~chromosome, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  labs(x = "Position on chromosome (Mbp)",
       y = expression(paste("Average ", italic(d[xy])))) +
  theme(panel.spacing.x = unit(0.01, "in")) +
  ggtitle("Shared dxy outliers 50kb windows")
write.csv(shared_outliers_50kb, "shared_dxy_outliers_50kb.csv", row.names = F)
png(filename = "shared_dxy_outliers_50kb.png",
    height = 2.5, width = 9, units = "in", res = 400)
b
dev.off()


#write bed file for shared outliers to ID genes in these regions
write_delim(shared_outliers_10kb %>%
              filter(pop_comparison == "dav_new") %>%
              select(c(chromosome, window_pos_1, window_pos_2, zscore)) %>%
              mutate(chromosome = gsub("scaf", "scaffold", chromosome)),
            "shared_dxy_outliers_10kb.bed",
            delim = '\t', col_names = F)


write_delim(shared_outliers_50kb %>%
              filter(pop_comparison == "dav_new") %>%
              select(c(chromosome, window_pos_1, window_pos_2, zscore)) %>%
              mutate(chromosome = gsub("scaf", "scaffold", chromosome)),
            "shared_dxy_outliers_50kb.bed",
            delim = '\t', col_names = F)





#### Determine if the number of shared outliers is more than expected
# Due to random chance

#Find number of outliers for each comparison
dn_outlier <- as.numeric(sum(outliers_10kb$pop_comparison == "dav_new"))
dr_outlier <- as.numeric(sum(outliers_10kb$pop_comparison == "dav_rup"))

#Find number of total windows for each comparison
dn_windows <- as.numeric(sum(filtered_10kb$pop_comparison == "dav_new"))
dr_windows <- as.numeric(sum(filtered_10kb$pop_comparison == "dav_rup"))

#Find number of shared outliers
count_shared_outliers <- nrow(shared_outliers_10kb)/2

# There are 118 for dn, and 166 for dr
# This is out of 32160 and 32140 windows, respectively
# There are 70 shared outliers
# In box 1, I pull 118 numbers from 32160
# In box 2, I pull 166 numbers from 32140
# What are the chances that 70 of these are the same number?

for(i in 1:100000){
  if(i == 1){
    countvec <- vector()
  }
  box1 <- sample.int(n = 32160, size = 118, replace = F)
  box2 <- sample.int(n = 32140, size = 166, replace = F)
  
  countvec <- c(countvec, sum(box1%in% box2))
}

t <- as.data.frame(countvec)
t$countvec <- as.factor(t$countvec)

ggplot(data = t, aes(x = countvec)) + 
  geom_histogram(stat = "count") +
  stat_count(position = position_stack(vjust = 0.25))

