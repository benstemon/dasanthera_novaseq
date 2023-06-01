library(tidyverse)
library(multcomp)
library(ggsignif)
library(data.table)
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_hybridzone_results/")

#find shared fdm outliers (among other things) between species that form hybrid zones
########################################################
#get list of plotting files in the directory
patternlist <- c("stats_10kb", "stats_5kb", "stats_1kb", "stats_500bp")

#find shared outliers within x bp of each other. Specify Z-score threshold
searchsize = 1000000
Zthresh = 3
saveplotfile = "5kb"

for (j in 1:length(patternlist)){
  filelist <- list.files(recursive = F, pattern = patternlist[j])
  namemod <- gsub("stats_", "", patternlist[j])
  
  #set up df 
  combdf <- data.frame()
  for (i in 1:length(filelist)){
    tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
      mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
    
    combdf <- rbind(combdf, tmpdata)
  }
  
  #filter out bad scafs
  #filter for tests of interest
  #flip sign on tests with car as P2 -- and switch names
  badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
  interest <- c('car_rup_dav','car_new_dav','car_new_rup')
  filterdf <- combdf %>%
    filter(!chr %in% badscafs) %>%
    filter(test %in% interest) %>%
    mutate(chr = gsub("fold", "", chr)) %>%
    mutate(midpoint = (windowStart + windowEnd) / 2)
  
  
  
  #TESTING AND PLOTTING
  #perform t-test to see if f_dm means differ significantly from 0
  test_dif0 <- filterdf %>%
    group_by(test) %>%
    summarize(
      t_statistic = t.test(f_dM, mu = 0, alternative = "two.sided")$statistic,
      p_value = t.test(f_dM, mu = 0, alternative = "two.sided")$p.value,
      mean_fdm = mean(f_dM),
      sd_fdm = sd(f_dM)
    )
  write.csv(test_dif0, paste(namemod, "fdm_vs_0.csv", sep="_"))
  
  #writing different style plot for fdm distribution
  sup1 <- ggplot(filterdf, aes(x = f_dM)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 0, col = "red", linetype = 2) +
    facet_grid(test~.) +
    theme_bw(base_size = 7.5) +
    ggtitle(bquote(paste("Distribution of ", italic(f[dM]), " -- ", .(namemod)))) +
    ylab("Window count") +
    xlab(expression(italic(f[dM])))
  
  ggsave(sup1, filename = paste("fdm_hz_distribution_",namemod,".png",sep=""),
         device = "png", width = 2, height = 3, units = "in", dpi = 400)
  
  
  #t-test to see if f_dM in CRD differs from CND
  ttest <- t.test(f_dM ~ test, data = filterdf,
                  subset = test %in% c("car_rup_dav", "car_new_dav"))
  sink(paste(namemod, "ttest_CRD vs CND.txt", sep="_"))
  print(ttest)
  sink()
  
  #tukey contrasts
  filterdf$group <- factor(filterdf$test)
  lm_model <- lm(f_dM ~ group, data = filterdf)
  tukey_results <- glht(lm_model, linfct = mcp(group = "Tukey"))
  summary(tukey_results)
  sink(paste(namemod, "tukey_test_CRD cs. CND.txt", sep="_"))
  print(summary(tukey_results))
  sink()
  
  
  
  #boxplot of fdM
  boxplot <- ggplot(filterdf, aes(x = test, y = f_dM)) +
    geom_violin() +
    geom_signif(comparisons =  combn(sort(unique(filterdf$test)), 2, simplify = F),
                step_increase = 0.1, test = "t.test", test.args = ) +
    stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                 colour = "red", width = 0.2) +
    labs(title = "Distribution of f_dM by test",
         x = "Test", y = "f_dM") +
    theme_bw()
  pdf(paste(namemod, "boxplot_fdm.pdf", sep="_"))
  print(boxplot)
  dev.off()
  
  
  ## IDENTIFY TOP FDM VALUES AND FIND SHARED "OUTLIERS"
  
  # Identify top 0.5% f_dM values for each test
  #CRD_fdm <- filterdf %>% filter(test == "car_rup_dav") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))
  #CND_fdm <- filterdf %>% filter(test == "car_new_dav") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))
  #CNR_fdm <- filterdf %>% filter(test == "car_new_rup") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))
  
  #Identify points >= Z-score threshold from the mean
  CRD_fdm <- filterdf %>% 
    filter(test == "car_rup_dav") %>% 
    mutate(zscore = ((f_dM - mean(f_dM))/sd(f_dM))) %>%
    filter(zscore >=Zthresh)
  CND_fdm <- filterdf %>% 
    filter(test == "car_new_dav") %>% 
    mutate(zscore = ((f_dM - mean(f_dM))/sd(f_dM))) %>%
    filter(zscore >=Zthresh)
  
  #SHARED OUTLIERS BETWEEN CRD AND CND
  shared_outliers <- data.frame()
  # loop over each row in CND_fdm
  for (i in 1:nrow(CND_fdm)) {
    # filter CRD_fdm to find any matching rows
    matched_rows <- CRD_fdm %>%
      filter(chr == CND_fdm$chr[i], abs(midpoint - CND_fdm$midpoint[i]) <= searchsize)
    
    # if there are any matching rows, rbind the hit and the match
    if (nrow(matched_rows) > 0) {
      shared_outliers <- rbind(shared_outliers, CND_fdm[i, ], matched_rows)
    }
    #filter duplicated rows
    shared_outliers <- unique(shared_outliers)
  }
  
  
  
  #PLOTTING
  if (length(shared_outliers) >0){
    #make plot of distribution of f_dm across chromosomes
    png(paste(namemod, "fdm-outliers.png", sep="_"),
        height = 3.5, width = 10, units = "in", res = 400)
    print(
      ggplot(filterdf, aes(x = midpoint/1000000, y = f_dM)) +
        geom_point(size = 0.2) +
        geom_point(data = shared_outliers, colour = "red", size = 0.5) +
        facet_grid(test~chr, scales = 'free_x', space = 'free_x')
    )
    dev.off()
    
    
    # more focused plot -- bird-bee
    focusplot <- filterdf %>%
      filter(chr %in% unique(shared_outliers$chr)) %>%
      filter(test %in% unique(shared_outliers$test))
    unique(focusplot$chr)
    
    #make rectangles for highlighting
    rects <- shared_outliers %>%
      arrange(chr,windowStart) %>%
      group_by(chr) %>%
      mutate(section = cumsum(c(TRUE, diff(midpoint) > 1e6))) %>%
      group_by(chr, section) %>%
      summarise(min = min(windowStart)/1000000, max = max(windowEnd)/1000000)
    
    write.csv(rects, file = paste(namemod, "shared_regions.csv", sep="_"))
    
    
    #make names for facet text
    groupnames <- c('car_new_dav'='new-dav', 'car_rup_dav'='rup-dav',
                    "scaf_1086"="scaf_1086" ,
                    "scaf_2687"="scaf_2687" ,
                    "scaf_1087"="scaf_1087" ,
                    "scaf_2532"="scaf_2532" ,
                    "scaf_2533"="scaf_2533" ,
                    "scaf_2686"="scaf_2686" ,
                    "scaf_2684"="scaf_2684" ,
                    "scaf_2685"="scaf_2685")
    
    if(namemod == saveplotfile){
      save(focusplot, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig5a_focusplot.ggplot")
      save(rects, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig5a_rects.ggplot")
      save(rects, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig5a_rects.ggplot")
      save(shared_outliers, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig5a_shared_outliers.ggplot")
    }
    
    
    png(paste(namemod, "shared_outliers-CRD-CND.png", sep="_"),
        width = length(unique(shared_outliers$chr))*1.5, height = 3, units = "in", res = 400)
    print(
      ggplot(focusplot, aes(x = midpoint/1000000, y = f_dM)) +
        geom_point(size = 0.2) +
        geom_point(data = shared_outliers, aes(x = midpoint/1000000, y = f_dM),
                   colour = "red", size = 0.2, alpha = 1) +
        geom_rect(data = rects, inherit.aes = FALSE,
                  aes(xmin = min, xmax = max,
                      ymin = -Inf, ymax = Inf),
                  fill = "blue", alpha = 0.2) +
        geom_hline(yintercept = 0, linetype = 2, colour = "red") +
        facet_grid(test~chr, scales = 'free_x', space = 'free_x',
                   labeller = as_labeller(groupnames)) +
        theme_bw() +
        xlab("Position on chromosome (Mbp)") +
        ylab(expression(italic(f[dM])))
    )
    dev.off()
    
    #finally, write these shared outliers to a final output file
    #these can be used to identify genes within the region of interest
    write.csv(shared_outliers, file = paste(namemod, 'shared_outliers.csv', sep="_"))
  }
  
}


########################################################





#PLOTTING FOR ALL THREE SHARED OUTLIERS
# more focused plot -- ALL SHARED OUTLIERS
focusplot_all <- filterdf %>%
  filter(chr %in% unique(shared_outliers_all$chr)) %>%
  filter(test %in% unique(shared_outliers_all$test))

#make rectangles for highlighting
rects_all <- shared_outliers_all %>%
  group_by(chr, test) %>%
  summarise(min = min(windowStart)/1000000, max = max(windowEnd)/1000000) %>%
  mutate(min = min(min)) %>%
  mutate(max = max(max))
write.csv(rects_all, file = "shared_regions_all.csv")


#make names for facet text
groupnames <- c('car_new_dav'='new-dav', 'car_rup_dav'='rup-dav', 'car_new_rup'='new-rup',
                "scaffold_1086"="scaffold_1086", "scaffold_2686"="scaffold_2686")


pdf("shared_outliers-all.pdf", width = 2.5, height = 4.2)
ggplot(focusplot_all, aes(x = midpoint/1000000, y = f_dM)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers_all, aes(x = midpoint/1000000, y = f_dM),
             colour = "red", size = 0.4, alpha = 1) +
  geom_rect(data = rects_all, inherit.aes = FALSE,
            aes(xmin = min, xmax = max,
                ymin = -Inf, ymax = Inf),
            fill = "blue", alpha = 0.2) +
  facet_grid(test~chr, scales = 'free_x', space = 'free_x',
             labeller = as_labeller(groupnames)) +
  theme_bw() +
  xlab("Position on chromosome (Mbp)") 
dev.off()











#
#
# COMPARING TO TESTS FROM ORIGINAL FDM
################################################################################
#
#
#
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_fullspecies_sliding_results/v1")

#get list of plotting files in the directory
filelist <- list.files(recursive = T, pattern = "10kb_10000_2500_plottingfile")

#set up df 
combdf <- data.frame()
for (i in 1:length(filelist)){
  tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
    mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
  
  combdf <- rbind(combdf, tmpdata)
}

#filter out bad scafs
#filter for tests of interest
#flip sign on tests with car as P2 -- and switch names
badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
interest <- c('car_rup_dav','new_car_dav','new_car_rup')
filterdf <- combdf %>%
  filter(!chr %in% badscafs) %>%
  filter(test %in% interest) %>%
  mutate(midpoint = (windowStart + windowEnd) / 2) %>%
  mutate(
    test = case_when(
      test == 'new_car_dav' ~ 'car_new_dav',
      test == 'new_car_rup' ~ 'car_new_rup',
      TRUE ~ test
    ),
    f_dM = if_else(
      test %in% c('car_new_dav', 'car_new_rup'),
      -1 * f_dM,
      f_dM
    )
  )



#TESTING AND PLOTTING
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_hybridzone_results/")

#perform t-test to see if f_dm means differ significantly from 0
test_dif0 <- filterdf %>%
  group_by(test) %>%
  summarize(
    t_statistic = t.test(f_dM, mu = 0, alternative = "two.sided")$statistic,
    p_value = t.test(f_dM, mu = 0, alternative = "two.sided")$p.value
  )
write.csv(test_dif0, "fdm_vs_0_ALLCARD.csv")


#t-test to see if f_dM in CRD differs from CND
ttest <- t.test(f_dM ~ test, data = filterdf,
                subset = test %in% c("car_rup_dav", "car_new_dav"))
sink("ttest_CRD vs CND ALLCARD.txt")
print(ttest)
sink()

#tukey contrasts
filterdf$group <- factor(filterdf$test)
lm_model <- lm(f_dM ~ group, data = filterdf)
tukey_results <- glht(lm_model, linfct = mcp(group = "Tukey"))
summary(tukey_results)
sink("tukey_test_CRD cs. CND ALLCARD.txt")
print(summary(tukey_results))
sink()

#boxplot of fdM
boxplot <- ggplot(filterdf, aes(x = test, y = f_dM)) +
  geom_violin() +
  geom_signif(comparisons =  combn(sort(unique(filterdf$test)), 2, simplify = F),
              step_increase = 0.1, test = "t.test", test.args = ) +
  stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
               colour = "red", width = 0.2) +
  labs(title = "Distribution of f_dM by test",
       x = "Test", y = "f_dM") +
  theme_bw()
pdf("boxplot_fdm_ALLCARD.pdf")
boxplot
dev.off()


## IDENTIFY TOP FDM VALUES AND FIND SHARED "OUTLIERS"

# Identify top 0.5% f_dM values for each test
CRD_fdm <- filterdf %>% filter(test == "car_rup_dav") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))
CND_fdm <- filterdf %>% filter(test == "car_new_dav") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))
CNR_fdm <- filterdf %>% filter(test == "car_new_rup") %>% arrange(desc(f_dM)) %>% slice_head(n = round(nrow(filterdf)*0.005))


#SHARED OUTLIERS BETWEEN CRD AND CND
#find shared outliers within x bp of each other
searchsize = 1000000
shared_outliers <- data.frame()
# loop over each row in CND_fdm
for (i in 1:nrow(CND_fdm)) {
  # filter CRD_fdm to find any matching rows
  matched_rows <- CRD_fdm %>%
    filter(chr == CND_fdm$chr[i], abs(midpoint - CND_fdm$midpoint[i]) <= searchsize)
  
  # if there are any matching rows, rbind the hit and the match
  if (nrow(matched_rows) > 0) {
    shared_outliers <- rbind(shared_outliers, CND_fdm[i, ], matched_rows)
  }
  #filter duplicated rows
  shared_outliers <- unique(shared_outliers)
}


# more focused plot -- bird-bee
focusplot <- filterdf %>%
  filter(chr %in% unique(shared_outliers$chr)) %>%
  filter(test %in% unique(shared_outliers$test))

#make rectangles for highlighting
rects <- shared_outliers %>%
  group_by(chr, test) %>%
  summarise(min = min(windowStart)/1000000, max = max(windowEnd)/1000000)
write.csv(rects, file = "shared_regions_ALLCARD.csv")

#make names for facet text
groupnames <- c('car_new_dav'='new-dav', 'car_rup_dav'='rup-dav',
                "scaffold_1086"="scaffold_1086", "scaffold_2686"="scaffold_2686")

pdf("shared_outliers-CRD-CND_ALLCARD.pdf", width = 4, height = 3)
ggplot(focusplot, aes(x = midpoint/1000000, y = f_dM)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "red") +
  geom_point(size = 0.2) +
  geom_point(data = shared_outliers, aes(x = midpoint/1000000, y = f_dM),
             colour = "red", size = 0.4, alpha = 1) +
  geom_rect(data = rects, inherit.aes = FALSE,
            aes(xmin = min, xmax = max,
                ymin = -Inf, ymax = Inf),
            fill = "blue", alpha = 0.2) +
  facet_grid(test~chr, scales = 'free_x', space = 'free_x',
             labeller = as_labeller(groupnames)) +
  theme_bw() +
  xlab("Position on chromosome (Mbp)") 
dev.off()

#finally, write these shared outliers to a final output file
#these can be used to identify genes within the region of interest
write.csv(shared_outliers, file = 'shared_outliers_ALLCARD.csv')


################################################################################





#
#
# USING FRUTICOSUS AS P1
################################################################################
#
#
#
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_fullspecies_sliding_results/v1")

#get list of plotting files in the directory
filelist <- list.files(recursive = T, pattern = "10kb_10000_2500_plottingfile")

#set up df 
combdf <- data.frame()
for (i in 1:length(filelist)){
  tmpdata <- read.table(filelist[i], header = T, as.is = T) %>%
    mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
  
  combdf <- rbind(combdf, tmpdata)
}

#filter out bad scafs
#filter for tests of interest
#flip sign on tests with car as P2 -- and switch names
badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
interest <- c('fru_dav_new','fru_dav_rup','new_rup_fru')
filterdf <- combdf %>%
  filter(!chr %in% badscafs) %>%
  filter(test %in% interest) %>%
  mutate(midpoint = (windowStart + windowEnd) / 2) %>%
  mutate(
    test = case_when(
      test == 'new_rup_fru' ~ 'fru_new_rup',
      TRUE ~ test
    ),
    f_dM = if_else(
      test %in% c('fru_new_rup'),
      -1 * f_dM,
      f_dM
    )
  )



#TESTING AND PLOTTING
#make plot of distribution of f_dm across chromosomes
ggplot(filterdf, aes(x = midpoint/1000000, y = f_dM)) +
  geom_point(size = 0.2) +
  facet_grid(test~chr, scales = 'free_x', space = 'free_x')



#perform t-test to see if f_dm means differ significantly from 0
test_dif0 <- filterdf %>%
  group_by(test) %>%
  summarize(
    t_statistic = t.test(f_dM, mu = 0, alternative = "two.sided")$statistic,
    p_value = t.test(f_dM, mu = 0, alternative = "two.sided")$p.value
  )


#t-test to see if f_dM in CRD differs from CND
t.test(f_dM ~ test, data = filterdf,
       subset = test %in% c("car_rup_dav", "car_new_dav"))


#tukey contrasts
filterdf$group <- factor(filterdf$test)
lm_model <- lm(f_dM ~ group, data = filterdf)
tukey_results <- glht(lm_model, linfct = mcp(group = "Tukey"))
summary(tukey_results)

adjusted_p_values <- summary(tukey_results, test = adjusted("bonferroni"))$test$pvalues
names(adjusted_p_values) <- levels(filterdf$group)


#boxplot of fdM
boxplot <- ggplot(filterdf, aes(x = test, y = f_dM)) +
  geom_violin() +
  geom_signif(comparisons =  combn(sort(unique(filterdf$test)), 2, simplify = F),
              step_increase = 0.1, test = "t.test", test.args = ) +
  labs(title = "Distribution of f_dM by test",
       x = "Test", y = "f_dM") +
  theme_bw()
boxplot

################################################################################



