setwd("~/Desktop/")
library(tidyverse)
library(cowplot)
library(readr)
library(utils)

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
toponumber = 3
tmp_weights <- read.delim("twisst_popspecific_v2/twisst_rup86_t2.weights.txt.gz",
                          sep = '\t', skip = toponumber)
windows_and_weights <- cbind(windows_and_weights, tmp_weights)


#now, reorder first by scaffold, then by start column
windows_and_weights <- arrange(windows_and_weights, scaffold, start)



#now these need to be written to new files for weights and windows, respectively
#first, windows:
newwindow <- windows_and_weights[,1:6]

#write windows file
#con <- gzfile("REORDERED_WINDOWS_10kbtrees.tsv.gz", "w")
#write.table(as.data.frame(newwindow), con,
#            sep = '\t', quote = F, row.names = F)
#close(con)


#next, the weights (and header for the weights file):
#starts at col 8 and continues to the end of the table
newweights <- windows_and_weights[,8:ncol(windows_and_weights)]
weightsheader <- as.data.frame(readLines("twisst_popspecific_v2/twisst_rup86_t2.weights.txt.gz",
                                         n = toponumber))
con <- gzfile("REORDERED_WEIGHTS_rup86_t2.tsv.gz", "w")
write.table(as.data.frame(weightsheader), con,
            quote = F, row.names = F, col.names = F)
write.table(as.data.frame(newweights), con,
            sep = '\t', quote = F, row.names = F)
close(con)
########################################


#setup for plotting with twisst functions:
#MAKE SURE FILES ARE IN CORRECT FORMAT! (SEE ABOVE)

#some notes:
#1. Remove badscafs
#2. Try to plot in ggplot rather than with functions given from source -- more flexibility
#3. See how smoothing of weightings is done and try to copy it for ggplotting
#4. Plot t1s and t2s together, with scaffolds split horizontally


#PLOTTING FOR TEST 1 VS TEST 2'S!
########################################
#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/results/twisst/source_plot_twisst.R")

#window data file
window_data_file <- "~/Desktop/twisst_popspecific_v2/REORDERED_WINDOWS_10kbtrees.tsv.gz"

#info for test 1 and test 2 weights
weights_file1 <- "~/Desktop/twisst_popspecific_v2/REORDERED_WEIGHTS_rup101_t1.tsv.gz"
weights_file2 <- "~/Desktop/twisst_popspecific_v2/REORDERED_WEIGHTS_rup101_t2.tsv.gz"


#generate twisst data frames
twisst_data1 <- import.twisst(weights_files=weights_file1,
                             window_data_files=window_data_file)

twisst_data2 <- import.twisst(weights_files=weights_file2,
                              window_data_files=window_data_file)


#generate list for good scaffolds
scaflist <- names(twisst_data1$lengths)
scaflist <- scaflist[!scaflist %in% 
                       c("scaf_2531", "scaf_2446", "scaf_2151", "scaf_1085")]


#set up new df to append to and eventually plot
newdf <- data.frame()

#set up test 1 df in the newdf
for (i in 1:length(scaflist)){
  twisst_subdata <- subset.twisst.by.regions(twisst_data1, scaflist[i])
  twisst_data_smooth <- smooth.twisst(twisst_subdata, span_bp = 2000000, spacing = 5000)
  
  #make new df
  tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
  tmpdf$chromosome <- scaflist[i]
  tmpdf$testno <- "test1"
  tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
  
  
  #bind to newdf
  newdf <- rbind(newdf, tmpdf)
}

#and now add the test 2 information so we can plot all together
for (i in 1:length(scaflist)){
  twisst_subdata <- subset.twisst.by.regions(twisst_data2, scaflist[i])
  twisst_data_smooth <- smooth.twisst(twisst_subdata, span_bp = 2000000, spacing = 5000)
  
  #make new df
  tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
  tmpdf$chromosome <- scaflist[i]
  tmpdf$testno <- "test2"
  tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
  
  
  #bind to newdf
  newdf <- rbind(newdf, tmpdf)
}


#generate plots of each topo and average weightings for the two tests
t1 <- subset.twisst.by.regions(twisst_data1, scaflist)
plot.twisst.summary(t1, cols = c("green", "orange", "blue"), lwd=3, cex=0.7)
a <- recordPlot()

t2 <- subset.twisst.by.regions(twisst_data2, scaflist)
plot.twisst.summary(t2, cols = c("green", "orange", "blue"), lwd=3, cex=0.7)
b <- recordPlot()


#plot output, facet grid by chromosome and test
c <- ggplot(newdf, aes(x = pos/1000000)) +
  geom_ribbon(aes(y = topo1, ymin = 0, ymax = topo1), col = "green", fill = "green", alpha = 0.25) + 
  geom_ribbon(aes(y = topo2, ymin = 0, ymax = topo2), col = "orange", fill = "orange", alpha = 0.25) + 
  geom_ribbon(aes(y = topo3, ymin = 0, ymax = topo3), col = "blue", fill = "blue", alpha = 0.25) + 
  ylim(c(0,1)) +
  facet_grid(testno ~ chromosome,
             scales = "free_x", space = "free_x") +
  xlab("Position on scaffold (Mb)") +
  ylab("Weightings")


#generate plot
rightpanel <- plot_grid(a, b, ncol = 1, scale = c(0.8, 0.8))

pdf("twisstplot_rup101_tests.pdf", width = 40, height = 10)
plot_grid(c, rightpanel, ncol = 2, rel_widths = c(6,1), scale = c(0.95, 0.95))
dev.off()

########################################





#PLOTTING DIFFERENCES ACROSS TAXA
########################################
setwd("~/Desktop/twisst_popspecific_v2/")

#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/source_plot_twisst.R")


#make list of comparison tests and basename for files to read in
basename = "REORDERED_WEIGHTS_"
testnames <- c("dav116", "new-dav", "rup86", "rup101")


#specify the window data file
window_data_file <- "~/Desktop/twisst_popspecific_v2/REORDERED_WINDOWS_10kbtrees.tsv.gz"


#make df to store all the relevant data
newdf <- data.frame()


#parameters to tune:
param_twisst_span = 2000000
param_twisst_space = 5000

#for loop which does the following:
#for each testname, read in the corresponding weights files
#for each of the main scaffolds, subset the data, and smooth at specified bp and spacing
#finally, add some additional ID information, and write to object "newdf"
for (j in 1:length(testnames)){
  #specify window and weights data files
  weights_file1 <- paste(basename, testnames[j], "_t1.tsv.gz", sep = "")
  weights_file2 <- paste(basename, testnames[j], "_t2.tsv.gz", sep = "")
  
  
  #generate twisst data frames
  twisst_data1 <- import.twisst(weights_files=weights_file1,
                                window_data_files=window_data_file)
  
  twisst_data2 <- import.twisst(weights_files=weights_file2,
                                window_data_files=window_data_file)
  
  #generate list for good scaffolds
  scaflist <- names(twisst_data1$lengths)
  scaflist <- scaflist[!scaflist %in% 
                         c("scaf_2531", "scaf_2446", "scaf_2151", "scaf_1085")]
  
  
  #set up test 1 df in the newdf
  for (i in 1:length(scaflist)){
    twisst_subdata <- subset.twisst.by.regions(twisst_data1, scaflist[i])
    twisst_data_smooth <- smooth.twisst(twisst_subdata,
                                        span_bp = param_twisst_span,
                                        spacing = param_twisst_space)
    
    #make new df, add extra information as needed
    tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
    tmpdf$chromosome <- scaflist[i]
    tmpdf$testno <- "test1"
    tmpdf$testname <- testnames[j]
    tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
    
    
    
    #bind to newdf
    newdf <- rbind(newdf, tmpdf)
  }
  
  #and now add the test 2 information so we can plot all together
  for (i in 1:length(scaflist)){
    twisst_subdata <- subset.twisst.by.regions(twisst_data2, scaflist[i])
    twisst_data_smooth <- smooth.twisst(twisst_subdata,
                                        span_bp = param_twisst_span,
                                        spacing = param_twisst_space)
    
    #make new df
    tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
    tmpdf$chromosome <- scaflist[i]
    tmpdf$testno <- "test2"
    tmpdf$testname <- testnames[j]
    tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
    
    
    #bind to newdf
    newdf <- rbind(newdf, tmpdf)
  }
}



#now we have a large df with all the plotting information.
#But we want to make an additional metric
#we want to know the difference between each of the topo values in test 1 vs test 2
#for each window, for each comparison
filtered_newdf <- newdf %>%
  group_by(pos, chromosome, testname) %>%
  mutate(diff_concordant = -diff(topo1)) %>%
  mutate(diff_disc1 = -diff(topo2)) %>%
  mutate(diff_disc2 = -diff(topo3)) %>%
  filter(!testno == "test2") %>%
  select(-c(testno, topo1, topo2, topo3)) %>%
  ungroup()

#restructure this to have an easier time plotting:
finalplot <- pivot_longer(data = filtered_newdf,
                     cols = c("diff_concordant", "diff_disc1", "diff_disc2"),
                     names_to = c("topology_diff"))


#now we can plot each of the four comparisons on a single plot
#to see how the differences in twisst topology weights may be similar across comparisons
legendcolors <- c("diff_concordant" = "green", "diff_disc1" = "orange", "diff_disc2" = "blue")

pdf("topodiffs_twisst_cross-comparisons.pdf", width = 15, height = 9)
ggplot(finalplot, aes(x = pos/1000000, y = value)) +
  geom_line(aes(col = topology_diff)) +
  facet_grid(testname ~ chromosome,
             scales = "free_x", space = "free_x") +
  xlab("Position on scaffold (Mb)") +
  ylab("Difference in weighting between test 1 and test 2") +
  scale_color_manual(values = legendcolors)
dev.off()
########################################
