setwd("~/Desktop/")
library(tidyverse)


#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/results/twisst/source_plot_twisst.R")

#for setting up window data file:
#windows as current are not ordered numerically, so there needs to be reordering
#both for the windows file and the weights files
########################################
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
toponumber = 15
tmp_weights <- read.delim("twisst_analyses/twisst_NCRD_10kbtrees.weights.txt.gz",
                          sep = '\t', skip = toponumber)
windows_and_weights <- cbind(windows_and_weights, tmp_weights)


#now, reorder first by scaffold, then by start column
windows_and_weights <- arrange(windows_and_weights, scaffold, start)



#now these need to be written to new files for weights and windows, respectively
library(readr)
library(utils)

#first, windows:
newwindow <- windows_and_weights[,1:6]

con <- gzfile("REORDERED_WINDOWS_10kbtrees.tsv.gz", "w")
write.table(as.data.frame(newwindow), con,
            sep = '\t', quote = F, row.names = F)
close(con)


#next, the weights (and header for the weights file):
#starts at col 8 and continues to the end of the table
newweights <- windows_and_weights[,8:ncol(windows_and_weights)]
weightsheader <- as.data.frame(readLines("twisst_analyses/twisst_NCRD_10kbtrees.weights.txt",
                                         n = toponumber))
con <- gzfile("REORDERED_WEIGHTS_NCRD_10kbtrees.tsv.gz", "w")
write.table(as.data.frame(weightsheader), con,
            quote = F, row.names = F, col.names = F)
write.table(as.data.frame(newweights), con,
            sep = '\t', quote = F, row.names = F)
close(con)
########################################


#setup for plotting with twisst functions:
#MAKE SURE FILES ARE IN CORRECT FORMAT! (SEE ABOVE)

#PLOTTING FOR NCR! (3 topologies)
########################################
weights_file <- "~/Desktop/twisst_analyses/REORDERED_WEIGHTS_NCR_10kbtrees.tsv.gz"
window_data_file <- "~/Desktop/twisst_analyses/REORDERED_WINDOWS_10kbtrees.tsv.gz"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


#get list of scaffold names
scaflist <- names(twisst_data$lengths)

setwd('~/Desktop/twisst_plotting_output_NCR/')

#loop to generate twisst plots for each scaffold
for (i in 1:length(scaflist)){
  #generate subset data
  twisst_subdata <- subset.twisst.by.regions(twisst_data, scaflist[i])
  
  #start pdf plot
  pdf(paste(scaflist[i], "twisst_plot_NCR.pdf", sep = '_'))
  
  #plot summary
  plot.twisst.summary(twisst_subdata, lwd=3, cex=0.7)
  
  #generate and plot smoothed data
  twisst_data_smooth <- smooth.twisst(twisst_subdata, span_bp = 2000000, spacing = 1000)
  plot.twisst(twisst_data_smooth, mode=2) #mode 2 overlays polygons, mode 3 would stack them
  plot.twisst(twisst_data_smooth, mode=3) #mode 2 overlays polygons, mode 3 would stack them
  dev.off()
}
########################################


#PLOTTING FOR NCRD! (15 topologies)
########################################
weights_file <- "~/Desktop/twisst_analyses/REORDERED_WEIGHTS_NCRD_10kbtrees.tsv.gz"
window_data_file <- "~/Desktop/twisst_analyses/REORDERED_WINDOWS_10kbtrees.tsv.gz"

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


#get list of scaffold names
scaflist <- names(twisst_data$lengths)


#keep only top 6 topologies
main_topologies <- order(twisst_data$weights_overall_mean, decreasing=T)[1:6]
twisst_data <- subset.twisst.by.topos(twisst_data, main_topologies)


setwd('~/Desktop/twisst_plotting_output_NCRD/')


#loop to generate twisst plots for each scaffold
for (i in 1:length(scaflist)){
  #generate subset data
  twisst_subdata <- subset.twisst.by.regions(twisst_data, scaflist[i])
  
  #start pdf plot
  pdf(paste(scaflist[i], "twisst_plot_NCRD.pdf", sep = '_'))
  
  #plot summary
  plot.twisst.summary(twisst_subdata, lwd=3, cex=0.7)
  
  #generate and plot smoothed data
  twisst_data_smooth <- smooth.twisst(twisst_subdata, span_bp = 2000000, spacing = 1000)
  plot.twisst(twisst_data_smooth, mode=2) #mode 2 overlays polygons, mode 3 would stack them
  plot.twisst(twisst_data_smooth, mode=3) #mode 2 overlays polygons, mode 3 would stack them
  dev.off()
}
########################################




