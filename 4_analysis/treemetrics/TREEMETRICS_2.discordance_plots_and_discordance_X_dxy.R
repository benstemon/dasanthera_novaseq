setwd('~/Desktop/')
library(ape)


#read in concatenated CDS trees
CDStrees <- read.tree("~/project storage/project_dasanthera_novaseq/results/trees/combined_CDStrees.tre")


#read in numbered CDS treepaths
treenumbers <- read.table("~/project storage/project_dasanthera_novaseq/results/treemetrics/numbered_CDStreepaths.txt")
colnames(treenumbers) <- c("number", "CDSname")
treenumbers[,2] <- 
  gsub(".fa.treefile", "", gsub("outfile_", "", treenumbers[,2]))


#read in bedfile with genome coordinates for CDS
#then match CDS names and paste genome coordinates
bedfile <- read.table("~/project storage/project_comparative_genome/genomes/annot_Pdavidsonii_1mb.gffread.genes.bed")

treenumbers$scaffold <- "NA"
treenumbers$bpstart <- "NA"
treenumbers$bpend <- "NA"

for (j in 1:nrow(bedfile)){
  if (bedfile[j,4] %in% treenumbers$CDSname == TRUE){
    treenumbers[match(bedfile[j,4], treenumbers$CDSname),][3:5] <- bedfile[j,1:3]
  }
}


#specify outgroup
out <- "CDS_bws_mon_61-7_S440.fixed.fa"
out2 <- "CDS_bws_lya_44-9_S439.fixed.fa"



#HERE is where customization comes into play
#####
#1. CDS_bws_rup_98-4_S426.fixed.fa -- CDS_bws_new_75-9_S431.fixed.fa
#2. CDS_bws_rup_98-4_S426.fixed.fa -- CDS_bws_new_80-5_S433.fixed.fa
#3. CDS_bws_rup_98-4_S426.fixed.fa -- CDS_bws_new_85-9_S432.fixed.fa
#4. CDS_bws_rup_86-8_S424.fixed.fa -- CDS_bws_new_75-9_S431.fixed.fa
#5. CDS_bws_rup_86-8_S424.fixed.fa -- CDS_bws_new_80-5_S433.fixed.fa
#6. CDS_bws_rup_86-8_S424.fixed.fa -- CDS_bws_new_85-9_S432.fixed.fa

#specify the individuals to analyze
car <- "CDS_bws_car_91-6_S436.fixed.fa"
rup <- "CDS_bws_rup_98-4_S426.fixed.fa"
new <- "CDS_bws_new_85-9_S432.fixed.fa"

#label new column for specific individuals included
#1. rup98_new75
#2. rup98_new80
#3. rup98_new85
#4. rup86_new75
#5. rup86_new80
#6. rup86_new85

newcolname <- "rup98_new85"
treenumbers[newcolname] <- "unknown"

#specify toy topologies, based on species tree topology
concordant <- read.tree(text = (paste("(",rup,",(",car,",",new,"));", sep = "")))
discordant1 <- read.tree(text = (paste("(",car,",(",rup,",",new,"));", sep = "")))
discordant2 <- read.tree(text = (paste("(",new,",(",rup,",",car,"));", sep = "")))
#####

#for each tree in the CDStrees file,
#1. root to the outgroup
#2. prune to the three focal tips
#3. test which tree toy tree is concordant
#4. log information
############################################################
for (i in 1:length(CDStrees)){
  #if outgroup taxon is in the tree, specify outgroup
  if(out %in% CDStrees[[i]]$tip.label | out2 %in% CDStrees[[i]]$tip.label){
    if(out %in% CDStrees[[i]]$tip.label){
      OG <- out
    }else{
      OG <- out2
    }
    
    #root tree with specified outgroup
    test_tree <- root(CDStrees[[i]], outgroup = OG)
    
    #if the specified tips are actually in the tree, prune it
    if(rup %in% test_tree$tip.label == TRUE &
       car %in% test_tree$tip.label == TRUE &
       new %in% test_tree$tip.label == TRUE){
      test_tree <- keep.tip(test_tree, c(rup,car,new))
      
      #if the pruned tree actually has branch lengths...
      if(sum(test_tree$edge.length) > 0){
        
        #compare pruned CDS tree topology to pruned species tree topology
        suppressWarnings(
          answer <- which(c(dist.topo(test_tree, concordant),
                            dist.topo(test_tree, discordant1),
                            dist.topo(test_tree, discordant2)) == 0))
        
        #determine CDS topology syntax and add to treenumbers df
        if(answer == 1){
          treenumbers[[newcolname]][i] <- ("concordant")
        }else if(answer == 2){
          treenumbers[[newcolname]][i] <- ("bird taxa sister")
        }else if(answer == 3){
          treenumbers[[newcolname]][i] <- ("discordant 2")
        }
      }else{
        #if no branch lengths, specify polytomy
        treenumbers[[newcolname]][i] <- ("polytomy")
      }
    }else{
      #if tips not in tree, specify missing ingroup 
      treenumbers[[newcolname]][i] <- "Missing ingroup"
    }
    
  }else{
    #if outgroup not in tree, specify missing outgroup
    treenumbers[[newcolname]][i] <- "Missing outgroup"
  }
}

#once finished with all relevant comparisons, write the file to a table
write.table(treenumbers, file = "NCR_CDS_topology_information.txt",
            quote = F, row.names = F, sep = "\t")
############################################################


##PART 2 -- generate discordance plots
############################################################
setwd("~/Desktop/")
library(tidyverse)
intable <- read.delim("NCR_CDS_topology_information.txt")


#"tidyversify" the input table
end_data <- data.frame()
for (i in 6:ncol(intable)){
  tmpdata <- intable[1:nrow(intable),c(1:5, i)]
  tmpdata$group <- colnames(tmpdata)[6]
  colnames(tmpdata)[6] <- "topology"
  end_data <- rbind(end_data, tmpdata)
}
remove(tmpdata)

#subset based on group
groupnames <- unique(end_data$group)

#for loop to produce plots for each unique group
pdf(file = "NCR_CDS_discordance_plots.pdf", width=10, height=12,pointsize=24)
for (i in 1:length(groupnames)){
  #generate subsetted data frame
  subdata <- end_data %>%
    filter(group == groupnames[i]) %>%
    filter(topology == "concordant" | 
             topology == "bird taxa sister" | 
             topology == "discordant 2") %>%
    mutate(bpstart = bpstart/1000000) %>%
    mutate(bpend = bpend/1000000)
  
  concord_prop <- sum(subdata$topology == "concordant")/ nrow(subdata)
  concord_prop <- round(concord_prop, digits = 3)
  
  bird_prop <- sum(subdata$topology == "bird taxa sister")/ nrow(subdata)
  bird_prop <- round(bird_prop, digits = 3)
  
  #plot
  print(
    ggplot(subdata, aes(x = bpstart, y = scaffold, colour = topology)) +
    geom_linerange(position = position_dodge(0.5), size = 2,
                   aes(xmin = bpstart, xmax = bpend, y = scaffold)) +
    xlab("Position on scaffold (Mb)") +
    ggtitle(paste(unique(subdata$group),
                  " concord prop = ", concord_prop, 
                  " bird prop = ", bird_prop,
                  sep = "")) +
    theme_bw()
  )
}
dev.off()
############################################################



##PART 3 -- generate discordance X Dxy plots
############################################################
setwd("~/Desktop/")
library(tidyverse)

#read in input table (topology information) and clean colnames a bit
intable <- read.delim("~/project storage/project_dasanthera_novaseq/results/treemetrics/NCR_CDS_topology_information.txt")
colnames(intable)[4:5] <- c("window_pos_1", "window_pos_2")


#"tidyversify" the input table
end_data <- data.frame()
for (i in 6:ncol(intable)){
  tmpdata <- intable[1:nrow(intable),c(1:5, i)]
  tmpdata$group <- colnames(tmpdata)[6]
  colnames(tmpdata)[6] <- "topology"
  end_data <- rbind(end_data, tmpdata)
}
remove(tmpdata)
remove(intable)



#read in dxy metrics for pop-specific comparisons
indxy <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixy_CDS/NCR_popspecific_dxy.txt")

#specify populations to include in the next steps
popset = c("rup_98", "rup_86",
           "new_75", "new_80", "new_85")

indxy <- indxy %>%
  filter(pop1 %in% popset) %>%
  filter(pop2 %in% popset) %>%
  mutate(group = paste(pop2, pop1, sep = "_")) %>%
  mutate(group = gsub("new_", "new", group)) %>%
  mutate(group = gsub("rup_", "rup", group)) %>%
  filter(!group %in% c("new80_new75",
                       "rup98_rup86",
                       "new80_new85",
                       "new75_new85"))

#ensure the filtering was successful -- groups are now identical across tables
unique(indxy$group) %in% unique(end_data$group)
unique(end_data$group) %in% unique(indxy$group)

#define all of the unique groups in the dfs
groupset <- unique(indxy$group)


#begin for loop -- a nested for loop
#first loop goes through each unique value in groupset to subset a data into group
for (i in 1:length(groupset)){
  #set up subsetted datasets for generating plots
  subdata_topology <- end_data %>%
    filter(group == groupset[i])
  subdata_dxy <- indxy %>%
    filter(group == groupset[i])
  
  #set up df for plotting
  plotting_product <- data.frame()
  
  #dirty nested loop to match windows and scaffold from topology df to dxy df
  #I am 100% sure there is a better way to do this, but it works
  for (j in 1:nrow(subdata_topology)){
    win1 <- subdata_topology$window_pos_1[j]
    win2 <- subdata_topology$window_pos_2[j]
    scaf <- subdata_topology$scaffold[j]
    
    
    tmpbind <- subdata_dxy[which(subdata_dxy$window_pos_1 == win1 &
                                   subdata_dxy$window_pos_2 == win2 &
                                   subdata_dxy$chromosome == scaf),]
    tmpbind <- cbind(tmpbind,
                     subdata_topology$number[j],
                     subdata_topology$CDSname[j],
                     subdata_topology$topology[j])
    
    plotting_product <- rbind(plotting_product, tmpbind)
  }
  
  #write the df generated from internal for loop
  colnames(plotting_product)[12:14] <- c("number", "CDSname", "topology")
  write.table(plotting_product,
              file = paste("merged_dxy_CDS-discordance_NCR_",
                           unique(plotting_product$group),
                           "_unfiltered.txt", sep = ""),
              quote = F, sep = "\t", row.names = F)
  
  #additional filters to exclude missing data, uninteresting topologies, etc.
  filtered_plotting_product  <- plotting_product %>%
    filter(!is.na(avg_dxy)) %>%
    filter(no_sites > 100) %>%
    #filter(count_missing/count_comparisons < 0.5) %>%
    filter(!topology %in% c("Missing ingroup", "Missing outgroup", "polytomy"))
  
  #plot the output
  pdf(file = paste("boxplots_dxy_vs_topologies_NCR_",
                   unique(filtered_plotting_product$group),
                   ".pdf",
                   sep = ""))
  print(ggplot(filtered_plotting_product, aes(x = topology, y = avg_dxy)) +
          #geom_violin() +
          geom_boxplot() +
          ggtitle(paste("dxy in different CDS topologies:",
                        unique(filtered_plotting_product$group)))
  )
  dev.off()
}
############################################################



#Part 4. Relative node depth metrics for the same CDS regions.
############################################################
setwd("~/Desktop/")
library(tidyverse)

#read in input table (topology information) and clean colnames a bit
intable <- read.delim("~/project storage/project_dasanthera_novaseq/results/treemetrics/NCR_CDS_topology_information.txt")
colnames(intable)[4:5] <- c("window_pos_1", "window_pos_2")

#"tidyversify" the input table
end_data <- data.frame()
for (i in 6:ncol(intable)){
  tmpdata <- intable[1:nrow(intable),c(1:5, i)]
  tmpdata$group <- colnames(tmpdata)[6]
  colnames(tmpdata)[6] <- "topology"
  end_data <- rbind(end_data, tmpdata)
}
remove(tmpdata)
remove(intable)

#read in dxy metrics for pop-specific comparisons
indxy <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixy_CDS/NCR_popspecific_dxy.txt")

#specify populations to include in the next steps
popset = c("rup_98", "rup_86",
           "new_75", "new_80", "new_85",
           "car_91")

indxy <- indxy %>%
  filter(pop1 %in% popset) %>%
  filter(pop2 %in% popset) %>%
  mutate(group = paste(pop2, pop1, sep = "_")) %>%
  mutate(group = gsub("new_", "new", group)) %>%
  mutate(group = gsub("rup_", "rup", group)) %>%
  mutate(group = gsub("car_", "car", group)) %>%
  filter(!group %in% c("new80_new75",
                       "rup98_rup86",
                       "new80_new85",
                       "new75_new85"))


#check the unique groups
#car_91 is involved with every topology, so group only reflects the new-car dichotomy
unique(indxy$group)



#when we do our filtering, we need three different groups in dxy:
#(1). new-car
#(2). new-rup
#(3). car-rup
#we then match each specific topology window to to all three values
#FOR loop to do all of this --- 

#define all of the unique groups in the dfs
groupset <- unique(indxy$group)

#we only want the first 6 groups in the groupset
for (i in 1:6){

  #specify which specific dxy comparisons you want for this trio of dxy
  tmp_groupset <- c(groupset[i],
                    paste(str_split(groupset[i], "_")[[1]][1], "car91", sep = "_"),
                    paste(str_split(groupset[i], "_")[[1]][2], "car91", sep = "_"))
  
  #now filter to only include these groupsets
  #set up subsetted datasets for generating plots
  #topology must be only first group
  subdata_topology <- end_data %>%
    filter(group == tmp_groupset[1])
  
  #dxy should have all three groups
  subdata_dxy <- indxy %>%
    filter(group %in% tmp_groupset)
  
  #dxy still needs a bit of filtering.
  #related rows are in sets of three. We want to keep 1 whole row, in addition to:
  #(1) remove the avg_dxy, count_missing, pop1, pop2, and group columns
  #(2) save and rename count_diffs and count_comparisons for each row (r_n, n_c, r_c)
  subdata_dxy <- subdata_dxy %>%
    select(-count_diffs, -pop1, -pop2, -group) %>%
    mutate(rn = ceiling(row_number()/3)) %>%
    group_by(chromosome, window_pos_1, window_pos_2, rn) %>%
    summarise(avg_dxy = paste0(avg_dxy, collapse = "," ),
              count_missing = paste0(count_missing, collapse = "," ),
              count_comparisons = paste0(count_comparisons, collapse = "," )) %>%
    ungroup() %>%
    arrange(rn) %>%
    separate(avg_dxy, c("dxy_r_n", "dxy_n_c", "dxy_r_c"), sep = ",") %>%
    separate(count_missing, c("count_missing_r_n", "count_missing_n_c", "count_missing_r_c"), sep = ",") %>%
    separate(count_comparisons, c("count_comparisons_r_n", "count_comparisons_n_c", "count_comparisons_r_c"), sep = ",")
  
  #set up df for plotting
  plotting_product <- data.frame()
  
  #dirty nested loop to match windows and scaffold from topology df to dxy df
  #I am 100% sure there is a better way to do this, but it works
  for (j in 1:nrow(subdata_topology)){
    win1 <- subdata_topology$window_pos_1[j]
    win2 <- subdata_topology$window_pos_2[j]
    scaf <- subdata_topology$scaffold[j]
    
    
    tmpbind <- subdata_dxy[which(subdata_dxy$window_pos_1 == win1 &
                                   subdata_dxy$window_pos_2 == win2 &
                                   subdata_dxy$chromosome == scaf),]
    tmpbind <- cbind(tmpbind,
                     subdata_topology$number[j],
                     subdata_topology$CDSname[j],
                     subdata_topology$topology[j])
    
    plotting_product <- rbind(plotting_product, tmpbind)
  }
  
  #write the df generated from internal for loop
  colnames(plotting_product)[14:16] <- c("number", "CDSname", "topology")
  
  #convert some columns from characters to numeric
  plotting_product[, 5:13] <- sapply(plotting_product[, 5:13], as.numeric)
  
  
  
  #some further filters on plotting data to remove uninteresting/not useful sites
  #different filters could be applied later. But I want to write at least this to df.
  #Also add calculation for RND, which is dependent on topology
  filtered_plotting_product  <- plotting_product %>%
    filter(!topology %in% c("Missing ingroup", "Missing outgroup", "polytomy")) %>%
    drop_na() %>%
    mutate(numerator = case_when(topology == "bird taxa sister" ~ dxy_r_n,
                                 topology == "concordant" ~ dxy_n_c,
                                 topology == "discordant 2" ~ dxy_r_c,)) %>%
    rowwise() %>%
    mutate(denominator = case_when(topology == "bird taxa sister" ~ mean(c(dxy_r_c, dxy_n_c)),
                                   topology == "concordant" ~ mean(c(dxy_r_c, dxy_r_n)),
                                   topology == "discordant 2" ~ mean(c(dxy_r_n, dxy_n_c)))) %>%
    mutate(RND = numerator/denominator)
  
  
  #write dxy and RND data to table
  write.table(filtered_plotting_product,
              file = paste("RNDPLOT_merged_dxy_CDS-discordance_NCR_",
                           unique(subdata_topology$group),
                           "_BASICFILTER.txt", sep = ""),
              quote = F, sep = "\t", row.names = F)
  
  
  #plot the output
  pdf(file = paste("violin_RND_vs_topologies_NCR_",
                   tmp_groupset[1],
                   ".pdf",
                   sep = ""))
  print(ggplot(filtered_plotting_product, aes(x = topology, y = RND)) +
          geom_violin() +
          #geom_boxplot() +
          ggtitle(paste("RND in different CDS topologies:",
                        unique(filtered_plotting_product$group)))
  )
  dev.off()
}
############################################################


#Part 4.1 Same as above but playing with plotting
############################################################
setwd("~/Desktop/RNDplots/")
library(tidyverse)
library(matrixStats)
infiles <- list.files(pattern = "BASICFILTER.txt")

for (i in 1:length(infiles)){
  #read in data
  indata <- read.table(infiles[i], sep = '\t', header = T)
  
  #filter data more stringently
  badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
  test <- indata %>%
    filter(window_pos_2 - window_pos_1 > 1000) %>%
    mutate(plotpos = rowMeans(select(., window_pos_2:window_pos_1))/1000000) %>%
    filter(!chromosome %in% badscafs) %>%
    filter_all(all_vars(!is.infinite(.))) %>%
    rowwise() %>%
    mutate(sd_counts = sd(c(count_comparisons_r_n, count_comparisons_n_c, count_comparisons_r_c))) %>%
    mutate(CV = sd_counts/mean(c(count_comparisons_r_n, count_comparisons_n_c, count_comparisons_r_c))) %>%
    filter(CV < 0.05)
  
  #plot of distribution of RND per scaffold
  pdf(paste("other_RNDplots",
            strsplit(strsplit(infiles[i], "_NCR_")[[1]][2], "_BASICFILTER")[[1]][1],
            ".pdf", sep = ""))
  print(
    ggplot(test, aes(x = RND, colour = topology)) +
      geom_density() +
      facet_grid(chromosome ~ .,
                 scales = "free_y", switch = "x", space = "free_x")
  )
  
  #plot of distribution of RND per topology
  print(
    ggplot(test, aes(x = RND, colour = topology)) +
      geom_density() +
      facet_grid(topology ~ .,
                 scales = "free_y", switch = "x", space = "free_x")
  )
  
  #histogram version of RND per topology
  print(
    ggplot(test, aes(x = RND, colour = topology)) +
      geom_histogram(binwidth = 0.01) +
      facet_grid(topology ~ .,
                 scales = "free_y", switch = "x", space = "free_x")
  )
  
  dev.off()
  
}

#plot of points along scaffolds
ggplot(test, aes(x = plotpos, y = RND)) +
  geom_point(size = 0.5, alpha = 0.5, aes(colour = topology)) +
  facet_grid(chromosome ~ .,)
############################################################

#Part 5. More RND stuff, but this time with topology defined by dxy itself
#Rather than from the IQtree gene tree estimates
############################################################
setwd("~/Desktop/RNDplots/")
library(tidyverse)
library(matrixStats)
library(gtable)
library(gridExtra)
infiles <- list.files(pattern = "BASICFILTER.txt")

for (i in 1:length(infiles)){
  #read in data
  indata <- read.table(infiles[i], sep = '\t', header = T)
  
  #specify bad scaffolds and filter data
  badscafs <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
  test <- indata %>%
    filter(window_pos_2 - window_pos_1 > 1000) %>%
    filter(!chromosome %in% badscafs) %>%
    mutate(plotpos = rowMeans(select(., window_pos_2:window_pos_1))/1000000) %>%
    rowwise() %>%
    mutate(dxymin = which.min(c(dxy_r_n, dxy_n_c, dxy_r_c))) %>%
    mutate(dxy_topo = case_when(dxymin == 1 ~ "bird taxa sister",
                                dxymin == 2 ~ "concordant",
                                dxymin == 3 ~ "discordant 2")) %>%
    mutate(dxy_numerator = case_when(dxy_topo == "bird taxa sister" ~ dxy_r_n,
                                     dxy_topo == "concordant" ~ dxy_n_c,
                                     dxy_topo == "discordant 2" ~ dxy_r_c,)) %>%
    mutate(dxy_denominator = case_when(dxy_topo == "bird taxa sister" ~ mean(c(dxy_r_c, dxy_n_c)),
                                       dxy_topo == "concordant" ~ mean(c(dxy_r_c, dxy_r_n)),
                                       dxy_topo == "discordant 2" ~ mean(c(dxy_r_n, dxy_n_c)))) %>%
    mutate(dxy_RND = dxy_numerator/dxy_denominator) %>%
    mutate(sd_counts = sd(c(count_comparisons_r_n, count_comparisons_n_c, count_comparisons_r_c))) %>%
    mutate(CV = sd_counts/mean(c(count_comparisons_r_n, count_comparisons_n_c, count_comparisons_r_c))) %>%
    filter(CV < 0.05) %>%
    filter_all(all_vars(!is.infinite(.)))
  
  #set up pdf for plotting
  pdf(paste("RNDplots_dxy_as_topology_",
            strsplit(strsplit(infiles[i], "_NCR_")[[1]][2], "_BASICFILTER")[[1]][1],
            ".pdf", sep = ""))
  
  #set up plotting substrate
  
  
  #plot of distribution of RND per topology
  p1 <- ggplot(test, aes(x = dxy_RND, colour = dxy_topo)) +
    geom_density() +
    facet_grid(dxy_topo ~ .,
               scales = "free_y", switch = "x", space = "free_x")
  
  #histogram version of RND per topology
  p2 <- ggplot(test, aes(x = dxy_RND, colour = dxy_topo)) +
    geom_histogram(binwidth = 0.02) +
    facet_grid(dxy_topo ~ .,
               scales = "free_y", switch = "x", space = "free_x") +
    theme(legend.position = "None")
  
  print(grid.arrange(p1, p2, nrow = 1))
  dev.off()
}



   
