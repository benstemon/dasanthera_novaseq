library(tidyverse)
library(cowplot)
setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_introtests_results/")

#Main comparisons -- 4 main tests across scaffolds
#########################################################
#we have 4 different window bins to consider as well. So we will make a big loop
#for the four bins
windowsizes <- c("10000", "5kb", "1kb", "500bp")

for (i in 1:length(windowsizes)){
  ws <- windowsizes[i]
  
  #read in dfs. First I am just interested in TWISST comparison for each of the four tests
  #dav-new-fru118
  
  t1 <- read.table(paste("dav-new",
                         list.files("dav-new", pattern = glob2rx(paste("fru118*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t1$comparison = "fru118-(dav-new)"
  
  #dav116-ncr-fru118
  t2 <- read.table(paste("dav116",
                         list.files("dav116", pattern = glob2rx(paste("fru118*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t2$comparison = "fru118-(dav116-ncr)"
  
  #rup86-new-car
  t3 <- read.table(paste("rup86",
                         list.files("rup86", pattern = glob2rx(paste("car*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t3$comparison = "car-(new-rup86)"
  
  #rup101-fru-dav
  t4 <- read.table(paste("rup101",
                         list.files("rup101", pattern = glob2rx(paste("dav*",ws,"*", sep = ""))),
                         sep = "/")
                   , as.is = T, header = T)
  t4$comparison = "dav-(fru-rup101)"
  
  
  #combine all of these plots into a single data frame, filter out badscafs
  #and generate a midpoint for plotting
  badscafs <- c("scaffold_2531","scaffold_2446","scaffold_2151","scaffold_1085")
  
  comdf <- rbind(t1, t2, t3, t4) %>%
    filter(!chr %in% badscafs) %>%
    rowwise() %>%
    mutate(midpoint = mean(c(windowStart, windowEnd))) %>%
    ungroup()
  
  
  #generate a plot -- f_dM across scaffolds
  a <- ggplot(comdf, aes(x = midpoint/1000000, y = f_dM)) +
    geom_point(size = 0.3, alpha = 0.7) +
    facet_grid(comparison ~ chr,
               scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    theme(axis.text = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 18),
          axis.title = element_text(size = 14)) +
    xlab("Position on scaffold (Mb)")
  
  
  b <- ggplot(comdf, aes(x = f_dM)) +
    geom_histogram(bins = 100) +
    facet_grid(comparison ~ .) +
    geom_segment(aes(x = 0 , y = 0, xend = 0, yend = Inf), color = "red", linetype = "dashed") +
    theme(strip.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    ggtitle("Genome-wide fdM counts")
  
  
  #write plots to output in combined plot
  pdf(paste("fdM_plots_maintests_",ws,".pdf", sep = ""), width = 20, height = 12)
  print(
  plot_grid(a, b, ncol = 2,
            rel_widths = c(6, 1))
  )
  dev.off()
  
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

#plot results
#generate a plot -- f_dM across scaffolds
a <- ggplot(comdf, aes(x = midpoint/1000000, y = f_dM)) +
  geom_point(size = 0.3, alpha = 0.7) +
  facet_grid(allopatric ~ chr,
             scales = "free_x", space = "free_x") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme(axis.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 18),
        axis.title = element_text(size = 14)) +
  xlab("Position on scaffold (Mb)")


b <- ggplot(comdf, aes(x = f_dM)) +
  geom_histogram(bins = 100) +
  facet_grid(allopatric ~ .) +
  geom_segment(aes(x = 0 , y = 0, xend = 0, yend = Inf), color = "red", linetype = "dashed") +
  theme(strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  ggtitle("Genome-wide fdM counts")


#write plots to output in combined plot
pdf(paste("fdM_plots_consistency_comparison_", folder, ".pdf", sep = ""),
    width = 20, height = length(unique(comdf$allopatric))*3)
plot_grid(a, b, ncol = 2,
          rel_widths = c(6, 1))
dev.off()
#########################################################

