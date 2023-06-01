library(tidyverse)
library(cowplot)
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_recenthybrids_results/")

#Main comparisons -- 3 main tests across scaffolds
#########################################################

#major regions highlighted by the TWISST results:
rects_dav116 <- data.frame(chromosome = c("scaf_1087", "scaf_2533", "scaf_2686"),
                           start = c(18.5, 0, 0),
                           end = c(25, 10, 9),
                           comparison = c("dav116","dav116","dav116"))
save(rects_dav116, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig4_dav116rects.ggplot")


rects_rup101 <- data.frame(chromosome = c("scaf_2533", "scaf_2687"),
                           start = c(0, 0),
                           end = c(23, 15),
                           comparison = c("rup101","rup101"))
save(rects_rup101, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig4_rup101rects.ggplot")

rects_rup86 <- data.frame(chromosome = c("scaf_1086", "scaf_1086",
                                         "scaf_1087","scaf_1087",
                                         "scaf_2532","scaf_2533",
                                         "scaf_2684", "scaf_2685",
                                         "scaf_2686", "scaf_2687"),
                          start = c(0, 45, 0, 30, 12, 30, 0, 0, 7, 20),
                          end = c(10, 55, 10, 45, 27, 40, 8, 10, 18, 50),
                          comparison = c("rup86","rup86","rup86","rup86","rup86",
                                       "rup86","rup86","rup86","rup86","rup86"))
save(rects_rup86, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig4_rup86rects.ggplot")




#we have 3 different window bins to consider as well. So we will make a big loop
#for the four bins
windowsizes <- c("10000", "5kb", "1kb", "500bp")
savewindow = "500bp"

for (i in 1:length(windowsizes)){
  ws <- windowsizes[i]
  
  #read in dfs. First I am just interested in TWISST comparison for each of the four tests
  
  #dav116-ncr-fru118
  t1 <- read.table(paste("dav116",
                         list.files("dav116", pattern = glob2rx(paste("fru118*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t1$comparison = "dav116"
  
  
  
  #rup101-fru-dav
  t3 <- read.table(paste("rup101",
                         list.files("rup101", pattern = glob2rx(paste("dav*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t3$comparison = "rup101"
  
  
  #rup86-new-car
  t4 <- read.table(paste("rup86",
                         list.files("rup86", pattern = glob2rx(paste("car*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t4$comparison = "rup86"
  
  
  #combine all of these plots into a single data frame, filter out badscafs
  #and generate a midpoint for plotting
  badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")
  
  comdf <- rbind(t1, t3, t4) %>%
    filter(!chr %in% badscafs) %>%
    rename(chromosome = chr) %>%
    mutate(chromosome = gsub("fold", "", chromosome)) %>%
    rowwise() %>%
    mutate(midpoint = mean(c(windowStart, windowEnd))) %>%
    ungroup()
  
  #save to figure plot compendium if conditions met
  if (ws == savewindow){
    save(comdf, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig4_fdms.ggplot")
  }
  
  #make separate df for just the fdm outliers
  outlierdf <- comdf %>%
    group_by(comparison) %>%
    slice_max(order_by = f_dM, prop = 0.05)
  
  
  #generate a minor plot: distribution of fdm
  sup1 <- ggplot(comdf, aes(x = f_dM)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 0, col = "red", linetype = 2) +
    facet_grid(comparison~.) +
    theme_bw(base_size = 7.5) +
    ggtitle(bquote(paste("Distribution of ", italic(f[dM]), " -- ", .(ws)))) +
    ylab("Window count") +
    xlab(expression(italic(f[dM])))
  
  ggsave(sup1, filename = paste("fdm_distribution_",ws,".png",sep=""),
         device = "png", width = 2, height = 3, units = "in", dpi = 400)
  
  #t-test to see if distributions differ from 0 (evidence of introgression)
  test_dif0 <- comdf %>%
    group_by(comparison) %>%
    summarize(
      t_statistic = t.test(f_dM, mu = 0, alternative = "two.sided")$statistic,
      p_value = t.test(f_dM, mu = 0, alternative = "two.sided")$p.value,
      mean_fdm = mean(f_dM),
      sd_fdm = sd(f_dM)
    )
  write.csv(test_dif0, paste("t-test_fdm_vs_0-", ws, ".csv", sep=""))
  
  
  #generate a plot -- f_dM across scaffolds
  a <- ggplot(comdf, aes(x = midpoint/1000000, y = f_dM)) +
    #geom_point(data = outlierdf, size = 0.3, col = "red") +
    facet_grid(comparison ~ chromosome,
               scales = "free_x", space = "free_x") +
    
    geom_rect(data = rects_dav116, inherit.aes = FALSE,
              aes(xmin = start, xmax = end,
                  ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.5) +
    geom_rect(data = rects_rup101, inherit.aes = FALSE,
              aes(xmin = start, xmax = end,
                  ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.5) +
    geom_rect(data = rects_rup86, inherit.aes = FALSE,
              aes(xmin = start, xmax = end,
                  ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.5) +
    geom_point(size = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", color = "blue", se = F, linewidth = 0.5) +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    xlab("Position on scaffold (Mb)") +
    ylab(expression(italic(f[dM])))
  
  
  b <- ggplot(comdf, aes(x = f_dM)) +
    geom_histogram(bins = 100) +
    facet_grid(comparison ~ .) +
    geom_segment(aes(x = 0 , y = 0, xend = 0, yend = Inf), color = "red", linetype = "dashed") +
    theme_bw() +
    ggtitle(expression(paste(italic(f[dM])," genome-wide"))) +
    xlab(expression(italic(f[dM]))) +
    theme(strip.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 8)) 
  
  
  #write plots to output in combined plot
  ggplot_alternative <- function(){
    plot_grid(a, b, ncol = 2,
              rel_widths = c(6, 1))
  }
  ggsave(paste("fdM_plots_maintests_",ws,".png", sep = ""),
         ggplot_alternative(),
         width = 10,
         height = 3.5,
         dpi = 400,
         bg = "transparent")
  
  
  #write outliers to file -- .bed format for post-processing
  write_delim(outlierdf,
              file = paste("fdm_outliers_",ws,".bed", sep = ""),
              delim = '\t')
  
  
  
}


#########################################################




#Plotting other comparisons -- check how consistent results are with dif. allopatric pops
#########################################################
#just change the folder name to the correct folder and everything else should be automatic
folder = "rup101"
filelist <- list.files(folder)


#make the combined data frame for this allopatric population
comdf <- data.frame()
for (i in 1:length(filelist)){
  tmptable <- read.table(paste(folder, filelist[i], sep = "/"), header = T, as.is = T)
  tmptable$allopatric = str_split(filelist[i], '_')[[1]][1]
  
  comdf <- rbind(comdf, tmptable)
}

#remove bad scafs, make midpoints for plotting
badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")
comdf <- comdf %>%
  filter(!chr %in% badscafs) %>%
  rowwise() %>%
  mutate(midpoint = mean(c(windowStart, windowEnd))) %>%
  ungroup()

#make outlier df
outlierdf <- comdf %>%
  group_by(allopatric) %>%
  slice_max(order_by = f_dM, prop = 0.005)

#plot results
#generate a plot -- f_dM across scaffolds
a <- ggplot(comdf, aes(x = midpoint/1000000, y = f_dM)) +
  geom_point(size = 0.3, alpha = 0.7) +
  geom_point(data = outlierdf, size = 0.3, col = "red") +
  facet_grid(allopatric ~ chr,
             scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size = 8)) +
  xlab("Position on scaffold (Mb)")


b <- ggplot(comdf, aes(x = f_dM)) +
  geom_histogram(bins = 100) +
  facet_grid(allopatric ~ .) +
  geom_segment(aes(x = 0 , y = 0, xend = 0, yend = Inf), color = "red", linetype = "dashed") +
  theme_bw() +
  ggtitle("fdM genome-wide") +
  theme(strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 8)) 


#write plots to output in combined plot
ggplot_alternative <- function(){
  plot_grid(a, b, ncol = 2,
            rel_widths = c(6, 1))
}
ggsave(paste("fdM_plots_consistency_comparison_",folder,".png", sep = ""),
       ggplot_alternative(),
       width = 10,
       height = length(unique(comdf$allopatric)) * 1.5,
       dpi = 300,
       bg = "transparent")
#########################################################

#with 4 main tests across scaffolds
#########################################################
#we have 4 different window bins to consider as well. So we will make a big loop
#for the four bins
windowsizes <- c("10000", "5kb", "1kb", "500bp")

for (i in 1:length(windowsizes)){
  ws <- windowsizes[i]
  
  #read in dfs. First I am just interested in TWISST comparison for each of the four tests
  
  #dav116-ncr-fru118
  t1 <- read.table(paste("dav116",
                         list.files("dav116", pattern = glob2rx(paste("fru118*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t1$comparison = "t1:fru118-(dav116-ncr)"
  
  
  #dav-new-fru118
  t2 <- read.table(paste("dav-new",
                         list.files("dav-new", pattern = glob2rx(paste("fru118*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t2$comparison = "t2:fru118-(dav-new)"
  
  
  #rup101-fru-dav
  t3 <- read.table(paste("rup101",
                         list.files("rup101", pattern = glob2rx(paste("dav*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t3$comparison = "t3:dav-(fru-rup101)"
  
  
  #rup86-new-car
  t4 <- read.table(paste("rup86",
                         list.files("rup86", pattern = glob2rx(paste("car*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t4$comparison = "t4:car-(new-rup86)"
  
  
  #combine all of these plots into a single data frame, filter out badscafs
  #and generate a midpoint for plotting
  badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")
  
  comdf <- rbind(t1, t2, t3, t4) %>%
    filter(!chr %in% badscafs) %>%
    rowwise() %>%
    mutate(midpoint = mean(c(windowStart, windowEnd))) %>%
    ungroup()
  
  #make separate df for just the fdm outliers
  outlierdf <- comdf %>%
    group_by(comparison) %>%
    slice_max(order_by = f_dM, prop = 0.005)
  
  #generate a plot -- f_dM across scaffolds
  a <- ggplot(comdf, aes(x = midpoint/1000000, y = f_dM)) +
    geom_point(size = 0.3, alpha = 0.7) +
    geom_point(data = outlierdf, size = 0.3, col = "red") +
    facet_grid(comparison ~ chr,
               scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    theme_bw() +
    theme(text = element_text(size = 8)) +
    xlab("Position on scaffold (Mb)")
  
  
  b <- ggplot(comdf, aes(x = f_dM)) +
    geom_histogram(bins = 100) +
    facet_grid(comparison ~ .) +
    geom_segment(aes(x = 0 , y = 0, xend = 0, yend = Inf), color = "red", linetype = "dashed") +
    theme_bw() +
    ggtitle("fdM genome-wide") +
    theme(strip.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 8)) 
  
  
  #write plots to output in combined plot
  ggplot_alternative <- function(){
    plot_grid(a, b, ncol = 2,
              rel_widths = c(6, 1))
  }
  ggsave(paste("fdM_plots_maintests_",ws,".png", sep = ""),
         ggplot_alternative(),
         width = 10,
         height = 6,
         dpi = 300,
         bg = "transparent")
  
  
  #write outliers to file -- .bed format for post-processing
  write_delim(outlierdf,
              file = paste("fdm_outliers_",ws,".bed", sep = ""),
              delim = '\t')
  
  
  
}


#########################################################