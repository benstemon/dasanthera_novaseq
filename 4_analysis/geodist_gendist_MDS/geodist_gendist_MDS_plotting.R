setwd('~/Desktop/')

#1. For plotting geodist vs. gendist
library(tidyverse)
library(geosphere)
####################################################
#read in data frames
raw_indxy <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixy_CDS/allpops_popspecific_dxy.txt")
inlatlon <- read.delim("dasanthera_coords.txt")

#change dav_118 to fru_118 in dxy and latlon
raw_indxy$pop1[raw_indxy$pop1 == "dav_118"] <- "fru_118"
raw_indxy$pop2[raw_indxy$pop2 == "dav_118"] <- "fru_118"
inlatlon$pop[inlatlon$pop == "dav_118"] <- "fru_118"


#join data frames
indxy <- raw_indxy %>% 
  filter(!pop1 %in% c('mon_61','lya_44')) %>%
  filter(!pop2 %in% c('mon_61','lya_44')) %>%
  left_join(inlatlon, by = c("pop1" = "pop")) %>%
  rename(lat1 = lat, lon1 = lon) %>%
  left_join(inlatlon, by = c("pop2" = "pop")) %>%
  rename(lat2 = lat, lon2 = lon) %>%
  mutate(popcomb = paste(pop1, pop2, sep = '_')) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(geo_dist = distGeo(c(lon1, lat1), c(lon2, lat2)))

  
#this is a test df to ensure the summarization appears to work correctly
#test <- indxy %>%
#  ungroup() %>%
#  select(!c(lon1, lat1, lon2, lat2, pop1, pop2, count_missing, no_sites)) %>%
#  filter(popcomb == "car_91_new_80")


#this actually summarizes avg dxy correctly across genome-wide for specific comparisons
#and then adds a tag for intra- or inter-specific comparison, and for pop-specific comparison
badscaf <- c('scaffold_2531','scaffold_2446','scaffold_2151','scaffold_1085')
indxy2 <- indxy %>%
  ungroup() %>%
  #filter(!chromosome %in% badscaf) %>%
  filter(!pop1 %in% c('fru_106','fru_118')) %>%
  filter(!pop2 %in% c('fru_106','fru_118')) %>%
  group_by(popcomb, pop1, pop2, geo_dist) %>%
  summarize(genomewide_dxy = sum(count_diffs)/sum(count_comparisons)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(comptype = case_when((str_split(popcomb, '_')[[1]][1] == str_split(popcomb, '_')[[1]][3]) == "TRUE" ~ "intraspecific",
                              (str_split(popcomb, '_')[[1]][1] == str_split(popcomb, '_')[[1]][3]) == "FALSE" ~ "interspecific")) %>%
  mutate(specpair = paste(sort(c(str_split(pop1, '_')[[1]][1], str_split(pop2, '_')[[1]][1])[1]),
                          sort(c(str_split(pop1, '_')[[1]][1], str_split(pop2, '_')[[1]][1])[2]),
                          sep = "_"))

#further iteratively filter to get sets of plots for each specific taxon (less cluttered)
library(data.table)
species_set <- c("car", "rup", "new", "dav")


#plot geo_dist vs. dxy
for (i in 1:length(species_set)){
  indxy3 <- indxy2[indxy2$specpair %like% species_set[i],]
  assign(paste("p", i, sep = ""),
         ggplot(indxy3, aes(x = geo_dist, y = genomewide_dxy, color = specpair)) +
           geom_point() + 
           geom_smooth(method = "lm", formula = y ~ x,
                       alpha = 0.5, linewidth = 0.5, se = F) +
           ggtitle(paste("species-level:", species_set[i], sep = " ")) +
           theme(axis.text = element_blank())
  )
}

library(gridExtra)
pdf("geo_dist_by_genomewide_dxy.pdf")
grid.arrange(p1, p2, p3, p4)
dev.off()



#for nls
ggplot(data = iris, aes(x = Sepal.Length,  y = Petal.Length, color = Species)) +
  geom_point() +
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1)))
####################################################



#2. For performing MDS and plotting those results.
setwd('~/Desktop/')
library(tidyverse)
library(reshape2)
library(ggrepel)

#read in data
raw_indxy <- read.delim("~/project storage/project_dasanthera_novaseq/results/pixy_CDS/allpops_popspecific_dxy.txt")

#specify unique pop combinations and recalculate genome-wide dxy
indxy <- raw_indxy %>% 
  mutate(popcomb = paste(pop1, pop2, sep = '_')) %>%
  drop_na() %>%
  group_by(popcomb, pop1, pop2) %>%
  summarize(genomewide_dxy = sum(count_diffs)/sum(count_comparisons)) %>%
  ungroup() %>%
  select(-popcomb)

#can use this to test if you don't trust the code
#test <- raw_indxy %>%
#  mutate(popcomb = paste(pop1, pop2, sep = '_')) %>%
#  filter(popcomb == "car_28_car_91") %>%
#  drop_na()
#sum(test$count_diffs)/sum(test$count_comparisons)

#there is a slight problem, in that the pairwise comparisons are not equal
#wrt who is listed first or second. To fix this you can create identical matrix
#but swap order of pop1 and pop2. Then all possible combinations exist.
indxy2 <- indxy[,c(2,1,3)]; colnames(indxy2)[1:2] = c('pop1', 'pop2')
indxy <- rbind(indxy, indxy2)

#use reshape2 package to transform this into a pairwise matrix
dxymat <- acast(indxy, pop1 ~ pop2, value.var = "genomewide_dxy")

#write to csv
write.csv(dxymat, file="pairwise_dxy_mat_allpops.csv", row.names = T)





#read in distance matrix, convert NA to 0
dxymat <- read.csv('pairwise_dxy_mat_allpops.csv', row.names = 1)
dxymat[is.na(dxymat)] <- 0

#perform MDS with cmdscale function
mds1 <- as.data.frame(cmdscale(dxymat, k =2))

#clean up
colnames(mds1) <- c("Score1", "Score2")
mds1$indiv <- rownames(mds1)
mds1$pop <- rownames(mds1)
mds1 <- mds1 %>%
  rowwise() %>%
  mutate(species = case_when(pop == "dav_118" ~ "fru",
                             pop != "dav_118" ~ str_split(pop, '_')[[1]][1]))



#plot
a <- ggplot(mds1, aes(x = Score1, y = Score2, col = species)) +
  geom_point() +
  scale_color_manual(values = c("purple", "skyblue","blue", "grey", "black", "red", "pink"))+
  geom_label_repel(aes(label = pop),
                   box.padding = 0.5,
                   point.padding = 0,
                   force = 25,
                   show.legend  = F) +
  ggtitle("genome-wide MDS")

pdf("MDS_genomewide-dxy_allpops.pdf")
print(a)
dev.off()

