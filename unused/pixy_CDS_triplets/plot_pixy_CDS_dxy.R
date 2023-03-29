setwd('~/Desktop/')
library(tidyverse)

#fullspecies
##############
rawdata <- read.delim("NCR_fullspecies_dxy.txt", header = T)

#make additional colum to specify the species interaction
rawdata <- rawdata %>%
  mutate(inter = paste(pop1, pop2, sep = " x "))


filterdata <- rawdata %>%
  na.omit() %>%
  subset(no_sites > 100) %>%
  subset(count_missing/count_comparisons < 0.5)
  
pdf("fullspecies_dxy.pdf")
ggplot(filterdata, aes(x = inter, y = avg_dxy, color = inter)) +
  geom_boxplot() +
  ggtitle("Full species Dxy") +
  ylim(0, 0.3) +
  theme(legend.position = "none")
dev.off()
##############


#popspecific
#SHASTA-rupicola
##############
rawdata <- read.delim("NCR_popspecific_dxy.txt", header = T)

#make additional colum to specify the species interaction
rawdata <- rawdata %>%
  mutate(inter = paste(pop1, pop2, sep = " x "))

#list the poulations we want in this set of comparisons
#shasta-rupicola: the shasta rupicola vs all newberryi pops
newpops <- c("new_84", "new_85", "new_75", "new_80")
ruppops <- c("rup_86")
carpops <- c("car_91")


#start for-loop here:

pdf("shasta-rupicola CDS dxy.pdf")
for (i in 1:length(newpops)){
  #define the popset (3 population comparison) to be used
  popset <- subset(
    unique(rawdata$inter),
    sapply(unique(rawdata$inter), grepl, pattern = ruppops) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops) | 
      sapply(unique(rawdata$inter), grepl, pattern = ruppops) &
      sapply(unique(rawdata$inter), grepl, pattern = newpops[i]) | 
      sapply(unique(rawdata$inter), grepl, pattern = newpops[i]) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops)
  )
  
  #filter the dataset to only include comparisons in the popset
  shasta_rup <- rawdata %>%
    na.omit() %>%
    subset(no_sites > 100) %>%
    subset(count_missing/count_comparisons < 0.5) %>%
    subset(inter == popset[1] | 
             inter == popset[2] | 
             inter == popset[3])
  
  #generate box plot
  print(
    ggplot(shasta_rup, aes(x = inter, y = avg_dxy, color = inter)) +
      geom_boxplot() +
      ggtitle(paste("Shasta-rupicola x", newpops[i])) +
      ylim(0, 0.3) +
      theme(legend.position = "none")
        )
}
dev.off()
##############
  

#popspecific
#SHASTA-newberryi
##############
rawdata <- read.delim("NCR_popspecific_dxy.txt", header = T)

#make additional colum to specify the species interaction
rawdata <- rawdata %>%
  mutate(inter = paste(pop1, pop2, sep = " x "))

#list the poulations we want in this set of comparisons
#shasta-rupicola: the shasta rupicola vs all newberryi pops
newpops <- c("new_85")
ruppops <- c("rup_86", "rup_101", "rup_105", "rup_98")
carpops <- c("car_91")


#start for-loop here:

pdf("shasta-newberryi CDS dxy.pdf")
for (i in 1:length(ruppops)){
  #define the popset (3 population comparison) to be used
  popset <- subset(
    unique(rawdata$inter),
    sapply(unique(rawdata$inter), grepl, pattern = ruppops[i]) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops) | 
      sapply(unique(rawdata$inter), grepl, pattern = ruppops[i]) &
      sapply(unique(rawdata$inter), grepl, pattern = newpops) | 
      sapply(unique(rawdata$inter), grepl, pattern = newpops) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops)
  )
  
  #filter the dataset to only include comparisons in the popset
  shasta_new <- rawdata %>%
    na.omit() %>%
    subset(no_sites > 100) %>%
    subset(count_missing/count_comparisons < 0.5) %>%
    subset(inter == popset[1] | 
             inter == popset[2] | 
             inter == popset[3])
  
  #generate box plot
  print(
    ggplot(shasta_new, aes(x = inter, y = avg_dxy, color = inter)) +
      geom_boxplot() +
      ggtitle(paste("Shasta-newberryi x", ruppops[i])) +
      ylim(0, 0.3) +
      theme(legend.position = "none")
  )
}
dev.off()
##############


#popspecific
#noSHASTA
##############
rawdata <- read.delim("NCR_popspecific_dxy.txt", header = T)

#make additional colum to specify the species interaction
rawdata <- rawdata %>%
  mutate(inter = paste(pop1, pop2, sep = " x "))

#list the populations we want in this set of comparisons
#shasta-rupicola: the shasta rupicola vs all newberryi pops
newpops <- c("new_75", "new_80")
ruppops <- c("rup_98")
carpops <- c("car_91")


#start for-loop here:

pdf("rattlesnake-rupicola CDS dxy.pdf")
for (i in 1:length(newpops)){
  #define the popset (3 population comparison) to be used
  popset <- subset(
    unique(rawdata$inter),
    sapply(unique(rawdata$inter), grepl, pattern = ruppops) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops) | 
      sapply(unique(rawdata$inter), grepl, pattern = ruppops) &
      sapply(unique(rawdata$inter), grepl, pattern = newpops[i]) | 
      sapply(unique(rawdata$inter), grepl, pattern = newpops[i]) &
      sapply(unique(rawdata$inter), grepl, pattern = carpops)
  )
  
  #filter the dataset to only include comparisons in the popset
  shasta_rup <- rawdata %>%
    na.omit() %>%
    subset(no_sites > 100) %>%
    subset(count_missing/count_comparisons < 0.5) %>%
    subset(inter == popset[1] | 
             inter == popset[2] | 
             inter == popset[3])
  
  #generate box plot
  print(
    ggplot(shasta_rup, aes(x = inter, y = avg_dxy, color = inter)) +
      geom_boxplot() +
      ggtitle(paste("rattlesnake-rupicola x", newpops[i])) +
      ylim(0, 0.3) +
      theme(legend.position = "none")
  )
}
dev.off()


