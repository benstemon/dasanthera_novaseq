setwd('~/Desktop/')
library(tidyverse)

#prepare list of regions containing genes of interest
#genelist <- read.csv("~/project storage/project_dasanthera_novaseq/dasanthera_genelist.csv")

#list of scaffolds to remove
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")

#function to transform pixy files to long format
pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}


#genic fraction files:
######
#10kb
genicfractionfile10kb <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_10kbwin_10kbslide.bed",
                                    header = F, sep = " ",
                                    col.names = c("chromosome", "window_pos_1", "window_pos_2", "genic_fraction"))
genicfractionfile10kb <- genicfractionfile10kb %>%
  mutate(window_pos_1 = window_pos_1+1) %>%
  mutate(mid = ceiling((window_pos_1 + window_pos_2)/2))

#50kb
genicfractionfile50kb <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_50kbwin_50kbslide.bed",
                                    header = F, sep = " ",
                                    col.names = c("chromosome", "window_pos_1", "window_pos_2", "genic_fraction"))
genicfractionfile <- genicfractionfile %>%
  mutate(window_pos_1 = window_pos_1+1)
######

# goal 1. Plotting pi, fst, dxy across parent pops.
########################
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_parentpops/")

pixydata_joined10kb <- pixy_to_long(list.files(pattern = "parentpops_10kb"))
pixydata_joined50kb <- pixy_to_long(list.files(pattern = "parentpops_50kb"))

#set up labeller
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)

#filter data frames to exclude bad scafs, 
pixyplotting10kb <- pixydata_joined10kb %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!is.na(value)) %>%
  mutate(chrom_color_group = case_when(as.numeric(factor(chromosome)) %% 2 != 0 ~ "even",
                                       TRUE ~ "odd" )) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst"))

pixyplotting50kb <- pixydata_joined50kb %>%
  filter(!chromosome %in% badscafs) %>%
  filter(!is.na(value)) %>%
  mutate(chrom_color_group = case_when(as.numeric(factor(chromosome)) %% 2 != 0 ~ "even",
                                       TRUE ~ "odd" )) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst"))


#plot these in ggplot
png("pixplot_davidsonii_50kb.png", height = 3.5, width = 9, units = "in", res = 400)
pixyplotting50kb %>%
  filter(pop1 == "davidsonii" | (pop1 == "rupicola" & pop2 == "davidsonii")) %>%
ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0) +
  facet_grid(statistic ~ chromosome,
             scales = "free", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value)) +
  ylab("Statistic Value")+
  ggtitle("davidsonii -- 50kb") +
  scale_color_manual(values = c("grey50", "black"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x.top = element_text(size = 6),
        legend.position ="none",
        panel.spacing = unit(0.01, "in"))
dev.off()
########################



#goal 2. compare pi, fst, dxy in windows with Ti > 0.7
####################################################################
# idea being -- maybe there is selection on the shared topology regions
# problem is that I don't have newberryi

#get TWISST functions
source("~/project storage/project_dasanthera_novaseq/source_plot_twisst.R")

#specify infiles
setwd("~/project storage/project_dasanthera_novaseq/results/twisst_fullspecies/")
window_data_file <- "REORDERED_WINDOWS_10kbtrees.tsv.gz"

weightfilelist <- c("REORDERED_WEIGHTS_fullspecies_car_rup_dav.txt.gz",
                    "REORDERED_WEIGHTS_fullspecies_new_car_dav.txt.gz")


#make combined list for the selected twisst windows
newdf <- data.frame()
for (i in 1:length(weightfilelist)){
  twisst_data <- import.twisst(weights_files = weightfilelist[i],
                               window_data_files = window_data_file)
  
  #make new df, add extra information as needed
  for(j in names(twisst_data$window_data)){
    tmpdf <- as.data.frame(twisst_data$window_data[j])
    tmpdf <- cbind(tmpdf, as.data.frame(twisst_data$weights[j]))
    colnames(tmpdf)[1:6] = colnames(twisst_data$window_data[[1]])
    colnames(tmpdf)[7:9] = colnames(twisst_data$weights[[1]])
    tmpdf$testname <- weightfilelist[i]
    
    #bind to newdf
    newdf <- rbind(newdf, tmpdf)
  }
}
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_parentpops/")

#topo 1 is the same, but topo 2 means different things because of the order
# for CRD, topo2 = Tn. for CND, topo2 = Ti.
tCRD <- newdf %>%
  filter(testname %in% "REORDERED_WEIGHTS_fullspecies_car_rup_dav.txt.gz") %>%
  select(-lnL) %>%
  filter(!scaffold %in% gsub("fold", "", badscafs)) %>%
  rename(Tc = topo1, Tn = topo2, Ti = topo3)

tCND <- newdf %>%
  filter(testname %in% "REORDERED_WEIGHTS_fullspecies_new_car_dav.txt.gz") %>%
  select(-lnL) %>%
  filter(!scaffold %in% gsub("fold", "", badscafs)) %>%
  rename(Tc = topo1, Ti = topo2, Tn = topo3)


#now the question is ----- is it better to do pixy here with my parentpops?
#it's just unclear to me how having these parent pops actually helps anything
#i could JUST do dav-rup...
newpixy10kb <- pixyplotting10kb %>%
  mutate(mid = ceiling((window_pos_1+window_pos_2)/2)) %>%
  mutate(chromosome = gsub("fold", "", chromosome)) %>%
  rename(scaffold = chromosome) %>%
  select(-window_pos_1, -window_pos_2, -chrom_color_group)


#inner join these by scaffold and mid, and then we can plot.
merged_plotting_df <-
  inner_join(newpixy10kb, tCRD, by = c("scaffold", "mid")) %>%
  #filter(Ti >= 0.6) %>%
  filter(statistic %in% "avg_wc_fst")

#plotting for pi
png("Ti_vs_pi_PARENTPOPS.png", height = 3.5, width = 4.5, units = "in", res = 400)
ggplot(merged_plotting_df, aes(x = Ti, y = value, group = pop1)) +
  geom_point(aes(colour = pop1), size = 0.5) +
  geom_smooth(method = "loess", se = F, aes(colour = pop1), span = 0.7) +
  theme_bw() +
  xlab("Weight of introgresion topology") +
  ylab(expression(pi))
dev.off()

#plotting for dxy, fst
png("Ti_vs_fst_PARENTPOPS.png", height = 3.5, width = 4.5, units = "in", res = 400)
ggplot(merged_plotting_df, aes(x = Ti, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess", se = F, span = 0.7) +
  theme_bw() +
  xlab("Weight of introgresion topology") +
  ylab(expression(F[ST])) +
  ggtitle(expression(paste("Ti weight vs. ", F[ST], " -- dav x rup")))
dev.off()
####################################################################



#3. compare dxy in parent pops vs. dxy in same species, but not near hybrid zones
# I am just going to re-load in everything to make it cleaner.

#read in pixy data -- ALLOPATRIC POPS
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_fullgenome_fullspecies_10kb/")
t1 <- read.table("fullgenome_fullspecies_10kb_dxy.txt", head = T) %>%
  filter(pop1 %in% c("P_davidsonii", "P_rupicola")) %>%
  filter(pop2 %in% c("P_davidsonii", "P_rupicola")) %>%
  mutate(pop1 = gsub("P_", "", pop1),
         pop2 = gsub("P_", "", pop2))

t2 <- read.table("fullgenome_fullspecies_10kb_fst.txt", head = T) %>%
  filter(pop1 %in% c("P_davidsonii", "P_rupicola")) %>%
  filter(pop2 %in% c("P_davidsonii", "P_rupicola")) %>%
  mutate(pop1 = gsub("P_", "", pop1),
         pop2 = gsub("P_", "", pop2))

t3 <- read.table("fullgenome_fullspecies_10kb_pi.txt", head = T) %>%
  pivot_wider(id_cols = c(chromosome, window_pos_1), id_expand = T,
              values_from = count_missing, names_from = pop) %>%
  dplyr::select(-P_cardwellii, -P_newberryi, -P_montanus, -P_lyallii, -P_fruticosus) %>%
  rename(missing_dav = P_davidsonii, missing_rup = P_rupicola)

#combine these into a single data set. Then, perform the missing data/quality filters
dxy_notsympatric <- inner_join(t1, t2) %>%
  inner_join(., t3) %>%
  filter(!chromosome %in% badscafs) %>%
  mutate(chromosome = gsub("fold", "", chromosome)) %>%
  na.omit() %>%
  #filters for high degree of missingness
  mutate(misrate_dav = missing_dav/count_comparisons,
         zscore_misrate_dav = (misrate_dav-mean(misrate_dav))/sd(misrate_dav),
         misrate_rup = missing_rup/count_comparisons,
         zscore_misrate_rup = (misrate_rup-mean(misrate_rup))/sd(misrate_rup),
         misrate = count_missing/count_comparisons,
         zscore_misrate = (misrate-mean(misrate))/sd(misrate)) %>%
  filter(abs(zscore_misrate_dav) <2,
         abs(zscore_misrate_rup) <2,
         abs(zscore_misrate) <2,) %>%
  #filters for large difference in missingness between populations
  mutate(missdiff = missing_dav/missing_rup) %>%
  filter(missdiff <2 & missdiff >0.5) %>%
  #filter for for high proportion of snps/sites
  mutate(snp_to_site = no_snps/no_sites,
         zscore_snp_to_site = (snp_to_site-mean(snp_to_site))/sd(snp_to_site)) %>%
  filter(abs(zscore_snp_to_site) <2) %>%
  #convert to long format and keep only dxy as value of interest
  pivot_longer(cols = avg_dxy) %>%
  dplyr::select(chromosome, window_pos_1, window_pos_2, value) %>%
  mutate(population = "not_sympatric") %>%
  mutate(mid = ceiling((window_pos_1+window_pos_2)/2)) #%>%
  #dplyr::select(-window_pos_1, -window_pos_2)
  

#read in pixy data -- SYMPATRIC POPS
setwd("~/project storage/project_dasanthera_novaseq/results/PIXY/pixyout_parentpops/")
t1 <- read.table("pixy_parentpops_10kb_dxy.txt", head = T)
t2 <- read.table("pixy_parentpops_10kb_fst.txt", head = T)
t3 <- read.table("pixy_parentpops_10kb_pi.txt", head = T) %>%
  pivot_wider(id_cols = c(chromosome, window_pos_1), id_expand = T,
              values_from = count_missing, names_from = pop) %>%
  rename(missing_dav = davidsonii, missing_rup = rupicola)

#combine these into a single data set. Then, perform the missing data/quality filters
dxy_sympatric <- inner_join(t1, t2) %>%
  inner_join(., t3) %>%
  filter(!chromosome %in% badscafs) %>%
  mutate(chromosome = gsub("fold", "", chromosome)) %>%
  na.omit() %>%
  #filters for high degree of missingness
  mutate(misrate_dav = missing_dav/count_comparisons,
         zscore_misrate_dav = (misrate_dav-mean(misrate_dav))/sd(misrate_dav),
         misrate_rup = missing_rup/count_comparisons,
         zscore_misrate_rup = (misrate_rup-mean(misrate_rup))/sd(misrate_rup),
         misrate = count_missing/count_comparisons,
         zscore_misrate = (misrate-mean(misrate))/sd(misrate)) %>%
  filter(abs(zscore_misrate_dav) <2,
         abs(zscore_misrate_rup) <2,
         abs(zscore_misrate) <2,) %>%
  #filters for large difference in missingness between populations
  mutate(missdiff = missing_dav/missing_rup) %>%
  filter(missdiff <2 & missdiff >0.5) %>%
  #filter for for high proportion of snps/sites
  mutate(snp_to_site = no_snps/no_sites,
         zscore_snp_to_site = (snp_to_site-mean(snp_to_site))/sd(snp_to_site)) %>%
  filter(abs(zscore_snp_to_site) <2) %>%
  #convert to long format and keep only dxy as value of interest
  pivot_longer(cols = avg_dxy) %>%
  dplyr::select(chromosome, window_pos_1, window_pos_2, value) %>%
  mutate(population = "sympatric") %>%
  mutate(mid = ceiling((window_pos_1+window_pos_2)/2)) #%>%
  #dplyr::select(-window_pos_1, -window_pos_2)


#COMBINE THE TWO FILTERED DATASETS
#inner join these dfs and generate a z score
merged_dxy <- inner_join(dxy_notsympatric, dxy_sympatric, by = c("chromosome", "mid")) %>%
  mutate(dxy_diff = value.y - value.x) %>%
  mutate(zscore = (dxy_diff - mean(dxy_diff))/sd(dxy_diff))
save(merged_dxy, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/filtered_merged_dxy-diff_10kb.obj")

#identify outliers
outlier_dxy <- merged_dxy %>%
  filter(zscore >= 3)
write.table(outlier_dxy, file = "outlier_dxy.txt", sep = '\t',
            row.names = F, quote = F)
save(outlier_dxy, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/outlier_dxy-diff_10kb.obj")

#write some statistics for the mean dxy values, whether they differ, 
#and same for dxy-diff
#x is not sympatric
dxystats <- data.frame(row.names = c("allopatric", "sympatric"))
dxystats$mean_dxy <- c(mean(merged_dxy$value.x), mean(merged_dxy$value.y))
dxystats$sd_dxy <- c(sd(merged_dxy$value.x), sd(merged_dxy$value.y))
dxystats$pval_ttest <- t.test(merged_dxy$value.x, merged_dxy$value.y)$p.value
dxystats$mean_dxy_diff <- mean(merged_dxy$dxy_diff)
dxystats$sd_dxy_diff <- sd(merged_dxy$dxy_diff)
dxystats$number_of_outliers <- nrow(outlier_dxy)
write.table(dxystats, file = "dxystats.txt", sep = '\t',
            row.names = T, col.names = T, quote = F)




png("test.png", units = "in", width = 9, height = 1.5, res = 400)
ggplot() +
  geom_point(data=merged_dxy, aes(x = mid/1000000, y = value.x), col = "red", size = 0.5, alpha = 0.5) +
  geom_point(data=merged_dxy, aes(x = mid/1000000, y = value.y), col = "blue", size = 0.5, alpha = 0.5) +
  geom_smooth(data = merged_dxy, method = "loess", se = F, aes(x = mid/1000000, y = value.x), col = "grey", linewidth = 0.75) +
  geom_smooth(data = merged_dxy, method = "loess", se = F, aes(x = mid/1000000, y = value.y), col = "black", linewidth = 0.75) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  theme_bw()
dev.off()

#Now... try to specifically find these outlier dxy_diff values
#that are also in similar regions as the fdm outliers...


#distribution of dxy for sympatric and non-sympatric
dxy_rbind <- rbind(dxy_notsympatric, dxy_sympatric)
a <- ggplot(dxy_rbind, aes(y = value, x = population)) +
  geom_boxplot() +
  ylab(expression(paste("mean ", d[xy]))) +
  theme_bw()


#overall distribution of dxydiff
b <- ggplot(merged_dxy, aes(x = dxy_diff)) +
  geom_histogram(bins = 150) +
  theme_bw() +
  xlab(expression(paste(d[xy-(sympatric)], " - ", d[xy-(not-sympatric)])))


#genomic distribution of dxy_diff and outliers
png("PLOT_all_dxy_diff_outliers.png", height = 1.5, width = 9, units = "in", res = 400)
ggplot(merged_dxy, aes(x = mid/1000000, y = dxy_diff)) +
  geom_point(size = 0.2) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  geom_point(data = outlier_dxy, col = "red", size = 0.25) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  ylab(expression(paste(d[xy], "_diff"))) +
  xlab("Position on chromosome (Mbp)")
dev.off()

#plot of SELECT chromosomes, with outlier windows from fdm analyses included:
rects <- data.frame(chromosome = c("scaf_1086", "scaf_1086", "scaf_2686"),
                    min = c(49.819556, 53.453444, 11.935612),
                    max = c(50.43394, 53.928787, 15.864582))

c <- merged_dxy %>%
  filter(chromosome %in% c("scaf_1086", "scaf_2686")) %>%
  ggplot(., aes(x = mid/1000000, y = dxy_diff)) +
  geom_rect(data = rects, inherit.aes = FALSE,
            aes(xmin = min, xmax = max,
                ymin = -Inf, ymax = Inf),
            fill = "blue", alpha = 0.2) +
  geom_point(size = 0.2) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  geom_point(data = outlier_dxy %>% filter(chromosome %in% c("scaf_1086", "scaf_2686")),
             col = "red", size = 0.25) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  ylab(expression(paste(d[xy], "_diff"))) +
  xlab("Position on chromosome (Mbp)")
png("PLOT_select_dxy_diff_outliers.png", height = 1.5, width = 3, units = "in", res = 400)
c
dev.off()


#close-up interesting regions

#2686: 11.935612	15.864582
outlier_region1 <- merged_dxy %>%
  mutate(mid = mid/1000000) %>%
  filter(chromosome %in% c("scaf_2686"),
         mid >= 11.935612-1.25 & mid <= 15.864582+1.25)
d <- ggplot(outlier_region1, aes(x = mid, y = dxy_diff)) +
  annotate("rect", xmin = 11.935612, xmax = 15.864582, ymin = -Inf, ymax = Inf,
            fill = "blue", alpha = 0.2) +
  geom_line() +
  geom_point(data =outlier_dxy %>%
               filter(chromosome%in%"scaf_2686"),
             aes(x = mid/1000000, y = dxy_diff), colour = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


#1086: 49.819556	50.43394
outlier_region2 <- merged_dxy %>%
  mutate(mid = mid/1000000) %>%
  filter(chromosome %in% c("scaf_1086"),
         mid >= 49.819556-0.75 & mid <= 50.43394+0.75)
e <- ggplot(outlier_region2, aes(x = mid, y = dxy_diff)) +
  annotate("rect", xmin = 49.819556, xmax = 50.43394, ymin = -Inf, ymax = Inf,
            fill = "blue", alpha = 0.2) +
  geom_line() +
  geom_point(data =outlier_dxy %>%
               filter(chromosome%in%"scaf_1086",
                      mid>=(49.819556-0.75)*1000000 & mid <=(50.43394+0.75)*1000000),
             aes(x = mid/1000000, y = dxy_diff), colour = "red") +
  theme_bw() +
  ylab(expression(paste(d[xy], "_diff"))) +
  theme(axis.title.x = element_blank())


#1086: 53.453444	53.928787
outlier_region3 <- merged_dxy %>%
  mutate(mid = mid/1000000) %>%
  filter(chromosome %in% c("scaf_1086"),
         mid >= 53.453444-0.25 & mid <= 53.928787+1.25)
f <- ggplot(outlier_region3, aes(x = mid, y = dxy_diff)) +
  annotate("rect", xmin = 53.453444, xmax = 53.928787, ymin = -Inf, ymax = Inf,
            fill = "blue", alpha = 0.2) +
  geom_line() +
  geom_point(data =outlier_dxy %>%
               filter(chromosome%in%"scaf_1086",
                      mid>=54*1000000 & mid <=55*1000000),
             aes(x = mid/1000000, y = dxy_diff), colour = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


#arrange these 
library(patchwork)
patchplot <- (a | b) /
  c /
  (e | f | d)

png("dxy_diff_plots.png", units = "in", width = 8, height = 6, res = 400)
patchplot
dev.off()

#we want to find the regions where dxydif is the highest
#the idea here is that regions of high dxy in sympatric parents
#is consistent with reinforcement
#question is -- are we that shocked? We know D and fd outliers tend to
#cluster in regions of low dxy. But it has to be that way...
#so it's "regular" dxy in parents, and just very low in non-parents?
#need experimental confirmation I suppose?

