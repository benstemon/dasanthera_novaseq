setwd("~/project storage/project_dasanthera_novaseq/results/twisst_fullspecies/")
library(tidyverse)
library(cowplot)
library(readr)
library(utils)

#for setting up window data file:
#windows as current are not ordered numerically, so there needs to be reordering
#both for the windows file and the weights files
########################################
infilelist <- list.files(pattern = "fullspecies")

for (i in 1:length(infilelist)){
  #set up initial windows file
  tmp_windows <- read.delim('~/project storage/project_dasanthera_novaseq/results/treemetrics/numbered_10kbtreepaths.txt', header = F, row.names = 1)
  colnames(tmp_windows) = "pathname"
  windows_and_weights <- tmp_windows %>%
    rowwise() %>%
    mutate(scaffold = str_split(pathname, "/")[[1]][[1]]) %>%
    mutate(start = as.numeric(str_split(str_split(pathname, "bp_")[[1]][2], "-")[[1]][1])) %>%
    mutate(end = as.numeric(str_split(str_split(pathname, "-")[[1]][2], ".fa")[[1]][1])) %>%
    mutate(mid = ceiling(mean(c(start, end)))) %>%
    mutate(sites = end-start+1) %>%
    mutate(lnL = NA) %>%
    ungroup() %>%
    mutate(treenumber = row_number()) %>%
    select(-pathname)
  
  
  #set up initial weights file and add to windows_and_weights
  #toponumber = the number of topologies at the beginning of the weights file
  toponumber = 3
  tmp_weights <- read.delim(infilelist[i],
                            sep = '\t', skip = toponumber)
  windows_and_weights <- cbind(windows_and_weights, tmp_weights)
  
  
  #now, reorder first by scaffold, then by start column
  windows_and_weights <- arrange(windows_and_weights, scaffold, start)
  
  
  
  #now these need to be written to new files for weights and windows, respectively
  #first, windows:
  #newwindow <- windows_and_weights[,1:6]
  
  #write windows file
  #con <- gzfile("REORDERED_WINDOWS_10kbtrees.tsv.gz", "w")
  #write.table(as.data.frame(newwindow), con,
  #            sep = '\t', quote = F, row.names = F)
  #close(con)
  
  
  #next, the weights (and header for the weights file):
  #starts at col 8 and continues to the end of the table
  newweights <- windows_and_weights[,8:ncol(windows_and_weights)]
  weightsheader <- as.data.frame(readLines(infilelist[i],
                                           n = toponumber))
  con <- gzfile(paste("REORDERED_WEIGHTS_", infilelist[i], sep = ""), "w")
  write.table(as.data.frame(weightsheader), con,
              quote = F, row.names = F, col.names = F)
  write.table(as.data.frame(newweights), con,
              sep = '\t', quote = F, row.names = F)
  close(con)
}


########################################


#setup for plotting with twisst functions:
#MAKE SURE FILES ARE IN CORRECT FORMAT! (SEE ABOVE)

#some notes:
# want two new metrics for twisst plots:
# 1. Discordance imbalance: weight(disc1) - weight(disc2)
# 2. Total Discordance: weight(disc1) + weight(disc2)


#PLOTTING DIFFERENCES ACROSS TAXA
########################################
#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/source_plot_twisst.R")

#specify the window data file
window_data_file <- "REORDERED_WINDOWS_10kbtrees.tsv.gz"

#make list of weightfiles
weightfilelist <- list.files(pattern = "REORDERED_WEIGHTS")

#parameters to tune:
param_twisst_span = 2000000
param_twisst_space = 5000



#make df to store all the relevant data
newdf <- data.frame()

#for loop which does the following:
#for each testname, read in the corresponding weights files
#for each of the main scaffolds, subset the data, and smooth at specified bp and spacing
#finally, add some additional ID information, and write to object "newdf"
for (i in 1:length(weightfilelist)){
  twisst_data <- import.twisst(weights_files = weightfilelist[i],
                               window_data_files = window_data_file)
  
  scaflist <- names(twisst_data$lengths)
  scaflist <- scaflist[!scaflist %in% 
                         c("scaf_2531", "scaf_2446", "scaf_2151", "scaf_1085")]
  
  for (j in 1:length(scaflist)){
    twisst_subdata <- subset.twisst.by.regions(twisst_data, scaflist[j])
    twisst_data_smooth <- smooth.twisst(twisst_subdata,
                                        span_bp = param_twisst_span,
                                        spacing = param_twisst_space)
    
    #make new df, add extra information as needed
    tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
    tmpdf$chromosome <- scaflist[j]
    tmpdf$testname <- weightfilelist[i]
    tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
    
    #bind to newdf
    newdf <- rbind(newdf, tmpdf)
  }
}

#now we have a large df with all the plotting information.
#But we want to add our additional metrics for each window, for each test
# 1. Discordance imbalance: weight(disc1) - weight(disc2)
# 2. Total Discordance: weight(disc1) + weight(disc2)
filtered_newdf <- newdf %>%
  group_by(pos, chromosome, testname) %>%
  mutate(discordance_imbalance = topo2-topo3) %>%
  mutate(total_discordance = topo2+topo3) %>%
  ungroup() %>%
  mutate(testname = str_remove(testname, "REORDERED_WEIGHTS_fullspecies_")) %>%
  mutate(testname = str_remove(testname, ".txt.gz"))

#add the genic fraction file
genicfractionfile <- read.delim("~/project storage/project_dasanthera_novaseq/results/miscellaneous/genicfraction_2mbwin_5kbslide.bed",
                                header = F, sep = " ",
                                col.names = c("chromosome", "pos", "end", "genic_fraction"))
genicfractionfile <- genicfractionfile %>%
  select(-end) %>%
  mutate(pos = pos+1) %>%
  mutate(chromosome = str_replace(chromosome, "scaffold", "scaf")) %>%
  filter(chromosome %in% scaflist)

#full join these together and pivot
finalplot <- full_join(filtered_newdf, genicfractionfile) %>%
  select(-c(topo1, topo2, topo3)) %>%
  drop_na() %>%
  pivot_longer(data = .,
               cols = c("discordance_imbalance", "total_discordance", "genic_fraction"),
               names_to = c("metric"))

#finalplot <- pivot_longer(data = finalplot,
#                          cols = c("discordance_imbalance", "total_discordance", "genic_fraction"),
#                          names_to = c("metric"))


pdf("twisst_fullspecies_10triplets_2mbwin_5kbslide.pdf", height = 10, width = 15)
ggplot(finalplot, aes(x = pos/1000000, y = abs(value))) +
  geom_line(aes(col = metric)) +
  facet_grid(testname ~ chromosome,
             scales = "free_x", space = "free_x") +
  xlab("Position on scaffold (Mb)")
dev.off()





