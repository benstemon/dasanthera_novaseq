library(tidyverse)

#call in the twisst plotting functions provided by devs
source("~/project storage/project_dasanthera_novaseq/source_plot_twisst.R")

#specify the window data file
window_data_file <- "REORDERED_WINDOWS_10kbtrees.tsv.gz"

#specify weight file
weightfile <- "REORDERED_WEIGHTS_fullspecies_new_car_rup.txt.gz"


#USING SMOOTHED WEIGHTS:
#parameters to tune:
param_twisst_span = 2000000
param_twisst_space = 5000

#make df to store all the relevant data
newdf <- data.frame()

#import twisst data
setwd("~/project storage/project_dasanthera_novaseq/results/twisst_fullspecies/")
twisst_data <- import.twisst(weights_files = weightfile,
                             window_data_files = window_data_file)

#set up list of scaffolds to keep
scaflist <- names(twisst_data$lengths)
scaflist <- scaflist[!scaflist %in% 
                       c("scaf_2531", "scaf_2446", "scaf_2151", "scaf_1085")]

#for loop which subsets to each of the scaffolds to keep
#then smooths twisst data and adds it to the final output df (newdf)
for (j in 1:length(scaflist)){
  twisst_subdata <- subset.twisst.by.regions(twisst_data, scaflist[j])
  twisst_data_smooth <- smooth.twisst(twisst_subdata,
                                      span_bp = param_twisst_span,
                                      spacing = param_twisst_space)
  
  #make new df, add extra information as needed
  tmpdf <- as.data.frame(twisst_data_smooth$pos); colnames(tmpdf)[1] <- "pos"
  tmpdf$chromosome <- scaflist[j]
  tmpdf$testname <- "NCR"
  tmpdf <- cbind(tmpdf, as.data.frame(twisst_data_smooth$weights))
  
  #bind to newdf
  newdf <- rbind(newdf, tmpdf)
}

#Topology 1 is the concordant topology...
#topo1 (mon,((new,car),rup));
#topo2 (mon,((new,rup),car));
#topo3 (mon,(new,(car,rup)));

# Convert the data frame to long format, which simplifies plotting along chromosomes
plotting_long <- newdf %>%
  pivot_longer(., cols = c(topo1, topo2, topo3),
               values_to = 'weight',
               names_to = 'topology')

#write this to the R-object-compendium for future plotting
NCR_twisst <- plotting_long
save(NCR_twisst, file = "~/project storage/project_dasanthera_novaseq/plot-making-compendium/twisst_NCR_weights.ggplot")


# Make plot of this
a <- ggplot(plotting_long, aes(x = pos/1000000, y = weight, group = topology)) +
  geom_line(aes(col = topology), linewidth = 0.5) +
  facet_grid(~chromosome, scales = "free_x", space = "free_x") +
  theme_bw() +
  xlab("Position on scaffold (Mbp)") + 
  ylab("Topology weight") +
  scale_color_manual(values = c("orange", "blue", "red"),
                     labels = c("concordant", "new-rup-sister", "car-rup-sister")) +
  ggtitle("Smoothed topology weights new-car-rup triplet") +
  labs(color = "Topology") +
  theme(panel.spacing.x = unit(0.01, "in"))

png("twisst_NCR.png", res = 400, units = "in", height = 2, width = 10)
a
dev.off()
