library(tidyverse)
library(multcomp)
library(ggsignif)
library(data.table)



#find shared fdm outliers (among other things) between species that form hybrid zones
########################################################
#get list of plotting files in the directory
patternlist <- c("stats_10kb", "stats_5kb", "stats_1kb", "stats_500bp")

dir.create("plots")


for (j in 1:length(patternlist)){
  setwd("~/project storage/project_dasanthera_novaseq/results/Dstats_new-rup/")
  filelist <- list.files(recursive = F, pattern = patternlist[j])
  namemod <- gsub("stats_", "", patternlist[j])
  
  #set up df 
  for (i in 1:length(filelist)){
    combdf <- read.table(filelist[i], header = T, as.is = T) %>%
      mutate(test = str_split(basename(filelist[i]), "_localFstats")[[1]][1])
  }
  
  #filter out bad scafs
  #filter for tests of interest
  #flip sign on tests with car as P2 -- and switch names
  badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
  interest <- c('car_rup_dav','car_new_dav','car_new_rup')
  filterdf <- combdf %>%
    filter(!chr %in% badscafs) %>%
    mutate(chr = gsub("fold", "", chr)) %>%
    mutate(midpoint = ((windowStart + windowEnd) / 2)/1000000)
  
  
  
  #TESTING AND PLOTTING
  #perform t-test to see if f_dm means differ significantly from 0
  setwd("plots")
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
  
  ggsave(sup1, filename = paste("fdm_new-rup_distribution_",namemod,".png",sep=""),
         device = "png", width = 2.5, height = 2.5, units = "in", dpi = 400)
  
  
    ## IDENTIFY TOP FDM VALUES AND PLOT ALONG GENOME
  #Identify points >= Z-score threshold from the mean
  outlier_fdm <- filterdf %>% 
    mutate(zscore = ((f_dM - mean(f_dM))/sd(f_dM))) %>%
    filter(zscore >=4)
  
  #PLOTTING
    #make plot of distribution of f_dm across chromosomes
    png(paste(namemod, "fdm-outliers.png", sep="_"),
        height = 1.5, width = 10, units = "in", res = 400)
    
    ggplot(filterdf, aes(x = midpoint, y = f_dM)) +
        geom_point(size = 0.2) +
        geom_point(data = outlier_fdm, colour = "red", size = 0.5) +
        facet_grid(~chr, scales = 'free_x', space = 'free_x') +
      theme_bw() +
      labs(x = "Position on scaffold (Mbp)", y = expression(italic(f[dM]))) +
      theme(panel.spacing.x = unit(0.01, "in"))
    dev.off()
    
    
    
    # more focused plot -- zoom in to f3p5ph
    # make rectangles for highlighting
    f3p5ph <- data.frame(chromosome = "scaf_1086",
                         start = 40.497697,
                         end = 40.899408)
    
    #plot the f3p5ph scaffold
    png(paste(namemod, "f3p5ph_scaffold.png", sep = "_"),
        width = 4, height = 3, units = "in", res = 400)
    print(filterdf %>% filter(chr %in% "scaf_1086") %>%
      ggplot(., aes(x = midpoint, y = f_dM)) +
      geom_rect(data = f3p5ph, inherit.aes = FALSE,
                aes(xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf),
                fill = "red", alpha = 0.5) +
      geom_point(size = 0.3, alpha = 0.7) +
      geom_hline(yintercept = mean(filterdf$f_dM), color = "red", linetype = "dashed") +
      geom_smooth(method = "loess", color = "blue", se = F, linewidth = 0.5) +
      theme_bw() +
      labs(x = "Position on scaffold (Mbp)", y = expression(italic(f[dM]))))
    dev.off()
    
    #plot a zoomed in version of this
    png(paste(namemod, "f3p5ph_scaffold_zoom.png", sep = "_"),
        width = 4, height = 3, units = "in", res = 400)
    print(
    filterdf %>% filter(chr %in% "scaf_1086",
                        midpoint >= 38 & midpoint <= 44) %>%
      ggplot(., aes(x = midpoint, y = f_dM)) +
      geom_hline(yintercept = mean(filterdf$f_dM), color = "red", linetype = "dashed") +
      geom_rect(data = f3p5ph, inherit.aes = FALSE,
                aes(xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf),
                fill = "red", alpha = 0.5) +
      geom_point(size = 0.3, alpha = 0.7) +
      theme_bw() +
      labs(x = "Position on scaffold (Mbp)", y = expression(italic(f[dM])))
    )
    dev.off()
  
}


#save filtered data.frame to permanent.obj for future plotting
save(filterdf, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/new-rup-fdm_500bp.ggplot")






