setwd("~/project storage/project_dasanthera_novaseq/results/twisst_IB-test/")
library(tidyverse)
library(cowplot)
library(readr)
library(utils)

#STEP 1. STOP!
#Ensure this step has been performed before proceeding. for setting up window data file:
#windows as current are not ordered numerically, so there needs to be reordering
#both for the windows file and the weights files
########################################
infilelist <- list.files(pattern = "IBL-test")

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
  toponumber = 105
  tmp_weights <- read.delim(infilelist[i],
                            sep = '\t', skip = toponumber)
  windows_and_weights <- cbind(windows_and_weights, tmp_weights)
  
  
  #now, reorder first by scaffold, then by start column
  windows_and_weights <- arrange(windows_and_weights, scaffold, start)
  
  
  
  #now these need to be written to new files for weights and windows, respectively
  #first, windows:
  newwindow <- windows_and_weights[,1:6]
  
  #write windows file
  con <- gzfile("REORDERED_WINDOWS_10kbtrees.tsv.gz", "w")
  write.table(as.data.frame(newwindow), con,
              sep = '\t', quote = F, row.names = F)
  close(con)
  
  
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


#STEP 2. Identify topologies consistent with internal branches in the species tree
##########################################
#we have three internal branches.
#we will search the edge matrix for MRCA for each group of taxa
#if clade exists, the number of times the edge appears is 2n (n = # of tips)
#1. rup-new-car -- x = 4
#2. dav-fru -- x = 2
#3. new-car -- x = 2

#need to eventually sum together tree weights for each of these three focal branches.
library(ape)
library(matrixStats)

#read in file with tree topologies specified by TWISST.
#these are at the top of the REORDERED_WEIGHTS file
intrees <- read.tree("IB-test_trees.tre")


#loop to make list of trees with the specified IBs
#could generalize this to turn into a function...
for (i in 1:length(intrees)){
  #pull one tree
  looptree <- intrees[[i]]
  
  #find MRCA branches
  IB1 <- getMRCA(looptree, c("new", "car", "rup"))
  IB2 <- getMRCA(looptree, c("new", "car"))
  IB3 <- getMRCA(looptree, c("dav", "fru"))
  
  
  #if this is the first iteration, make the list to store information
  if(i == 1){
    IB_topolist <- list("IB1" = vector(), "IB2" = vector(), "IB3" = vector())
  }
  
  #if IB occurs in matrix specified number of times,
  #add the tree number to the list for that internal branch
  #IB1
  if(sum(rowCounts(mrca(looptree), value = IB1)) == 4){
    IB_topolist$IB1 <- append(IB_topolist$IB1, i)
  }
  
  #IB2
  if(sum(rowCounts(mrca(looptree), value = IB2)) == 2){
    IB_topolist$IB2 <- append(IB_topolist$IB2, i)
  }
  
  #IB3
  if(sum(rowCounts(mrca(looptree), value = IB3)) == 2){
    IB_topolist$IB3 <- append(IB_topolist$IB3, i)
  }
  
}

#now we have the list of topologies consistent with each of the focal internal
#branches in the species tree.
##########################################

#STEP 3. Plotting weights for topologies consistent with internal branches
##########################################
#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/source_plot_twisst.R")

#specify the window data file
window_data_file <- "REORDERED_WINDOWS_10kbtrees.tsv.gz"

#make list of weightfiles
weightfilelist <- list.files(pattern = "REORDERED_WEIGHTS")


#read in twisst data
twisst_data <- import.twisst(weights_files = weightfilelist[1],
                             window_data_files = window_data_file)

#make list of scaffolds to be used for visualization
scaflist <- names(twisst_data$lengths)
scaflist <- scaflist[!scaflist %in% 
                       c("scaf_2531", "scaf_2446", "scaf_2151", "scaf_1085")]
#IB1 = new,car,rup
#IB2 = new,car
#IB3 = dav.fru
########################################



#METHOD a: my own code. Gets me the correct values, but smoothing isn't as good?
########################################
#for loop to generate data frames
for(i in 1:length(scaflist)){
  
  #if the first iteration...
  if(i ==1){
    
    #make new df
    t <- as.data.frame(twisst_data$pos[scaflist[i]]); colnames(t) = "midpoint"
    
    #fill in the scaffold column
    t$scaffold <- scaflist[i]
    
    #fill in total concordance column, based on summed weights
    #of topologies identified from IB search
    t$concordance_IB1 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB1])
    t$concordance_IB2 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB2])
    t$concordance_IB3 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB3])
    
  }else{
    
    #otherwise, make a temporary data frame with same information for this loop
    tmp <- as.data.frame(twisst_data$pos[scaflist[i]]); colnames(tmp) = "midpoint"
    tmp$scaffold <- scaflist[i]
    tmp$concordance_IB1 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB1])
    tmp$concordance_IB2 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB2])
    tmp$concordance_IB3 <- rowSums((twisst_data$weights[[scaflist[i]]])[IB_topolist$IB3])
    
    #merge the tmp data frame and the original
    t <- rbind(t, tmp)
  }
}

#convert to long format
plotting_long <- t %>%
  pivot_longer(., cols = c(concordance_IB1, concordance_IB2, concordance_IB3),
               values_to = "concordance_proportion") %>%
  na.omit() %>%
  mutate(discordance_proportion = 1-concordance_proportion)

#write this to the R-object-compendium for future plotting
save(plotting_long, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig3a.twisst_IB_weights_methodA.ggplot")


#generate plot
a <- ggplot(plotting_long, aes(x = midpoint/1000000, y = concordance_proportion, group = name)) +
  geom_smooth(aes(y = concordance_proportion, col = name), se = F, method = "loess", span = 0.05, linewidth = 0.5) +
  facet_grid(~scaffold,
             scales = "free_x",
             space = "free_x") +
  theme_bw() +
  #xlim(0,1) +
  xlab("Position on scaffold (Mb)") +
  ylab("Topology weight") +
  scale_color_manual(values = c("blue", "orange", "darkgreen"),
                     labels = c("IB1", "IB2", "IB3")) +
  ggtitle("Smoothed topology weights for trees concordant with internal branches") +
  labs(color = "Internal branch")
png("twisst_IB_method-a_0.05span.png", units = "in", height = 2, width = 10, res =400)
a
dev.off()


########################################


#METHOD b: twisst code for combining topologies...
#here I am summing the smoothed weights for all topologies consistent with internal branches
########################################
for (i in 1:length(scaflist)){
  
  if (i == 1){
    newdf <- data.frame()
  }
  
  twisst_subdata <- subset.twisst.by.regions(twisst_data, scaflist[i])
  twisst_data_smooth <- smooth.twisst(twisst_subdata,
                                      span_bp = 2000000,
                                      spacing = 5000)
  
  #make new df, add extra information as needed
  tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
  tmpdf$chromosome <- scaflist[i]
  
  #make object for summed weights of topologies consistent with IBs, for each IB
  conc_IB1 <- rowSums(twisst_data_smooth$weights[[1]][IB_topolist$IB1])
  conc_IB2 <- rowSums(twisst_data_smooth$weights[[1]][IB_topolist$IB2])
  conc_IB3 <- rowSums(twisst_data_smooth$weights[[1]][IB_topolist$IB3])
  
  #bind these to the tmpdf
  tmpdf <- cbind(tmpdf, conc_IB1, conc_IB2, conc_IB3)
  
  #bind to newdf
  newdf <- rbind(newdf, tmpdf)
}

#pivot to long format
plotting_long <- newdf %>%
  pivot_longer(., cols = c(conc_IB1, conc_IB2, conc_IB3),
               values_to = 'concordance_weight',
               names_to = 'internal_branch')


#write this to the R-object-compendium for future plotting
save(plotting_long, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/twisst_IB_weights_methodB.ggplot")


#plot the newdf
png("twisst_IB_method-b_2MB_5kb.png", units = "in", height = 2, width = 10, res =400)
ggplot(plotting_long, aes(x = pos/1000000, y = concordance_weight, group = internal_branch)) +
  geom_line(aes(col = internal_branch), linewidth = 0.5) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  theme_bw() +
  xlab("Position on scaffold (Mb)") + 
  ylab("Topology weight") +
  scale_color_manual(values = c("blue", "orange", "darkgreen"),
                     labels = c("IB1", "IB2", "IB3")) +
  ggtitle("Smoothed topology weights for trees concordant with internal branches") +
  labs(color = "Internal branch")
dev.off()



#####################


#STEP 4: use the site concordance factors as a metric for internal branch SUPPORT
#####################
#the internal branches we are interested in are:
#NCR (IB1) = 30
#NC (IB2) = 26
#DF (IB3) = 32

#read in table, filter to relevant nodes, and rename accordingly
scf_tab <- read.table('~/project storage/project_dasanthera_novaseq/results/treemetrics/concordance_factors/10kbsource_10kbref.cf.stat_loci',
                      header = T) %>%
  filter(ID %in% c(30,26,32)) %>%
  rename(internal_branch = ID) %>%
  mutate(internal_branch = gsub("30", "IB1", internal_branch)) %>%
  mutate(internal_branch = gsub("26", "IB2", internal_branch)) %>%
  mutate(internal_branch = gsub("32", "IB3", internal_branch))


#also read in information for the tree names and informativeness of trees
#tidy these up so they can be joined with twisst data later on
treenames <- read.delim('~/project storage/project_dasanthera_novaseq/results/treemetrics/concordance_factors/10kbsource_cf_site_ids.txt',
                        header = T, sep = '\t') %>%
  select(c(-Type, -Seqs, -Sites, -Model)) %>%
  rename(PartID = Subset) %>%
  mutate(Name = gsub("scaffold", "scaf", Name)) %>%
  mutate(chromosome = str_extract(Name, "scaf_\\d+")) %>%
  filter(chromosome %in% scaflist) %>%
  mutate(pos = ceiling((as.numeric(str_extract(Name, "(?<=bp_)\\d+")) + 
                          as.numeric(str_extract(Name, "(?<=-)\\d+")))/2))



#join these two dfs by PartID, and change names of the 
scfdata_long <- full_join(scf_tab, treenames, by = "PartID") %>%
  na.omit() %>%
  pivot_longer(., cols = c(sC, sD1, sD2))



#plot these data
pdf("SCFs_IB_unsmoothed.pdf", height = 9, width = 15)
ggplot(scfdata_long, aes(x = pos/1000000, y = value, group = name)) +
  #geom_smooth(aes(col = name), method = "loess", span = 0.05, linewidth = 0.75) +
  geom_line(aes(col = name), linewidth = 0.75) +
  facet_grid(internal_branch~chromosome, scales = "free_x", space = "free_x") +
  theme_bw() +
  xlab("Position on scaffold (Mb)")
dev.off()


#second plot -- just concordant triplet. Use this one for figure.
scfdata_short <- scfdata_long %>%
  pivot_wider(values_from = value, names_from = name)
#save to object
save(scfdata_short, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig3b.ggplot")
#load("~/project storage/project_dasanthera_novaseq/plot-making-compendium/fig3b.ggplot")


b <- ggplot(scfdata_short, aes(x = pos/1000000, y = sC, group = internal_branch)) +
  geom_line(aes(col = internal_branch), linewidth = 0.5) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  theme_bw() +
  xlab("Position on scaffold (Mb)") +
  ylab("sC") +
  scale_color_manual(values = c("blue", "orange","darkgreen"),
                     labels = c("IB1", "IB2", "IB3")) +
  ggtitle("Site concordance factors for the concordant triplet") +
  labs(color = "Internal branch")

png("SCFs_SConly_IB_0.05span.png", units = "in", height = 2, width = 10, res =400)
b
dev.off()

#####################
  
         

