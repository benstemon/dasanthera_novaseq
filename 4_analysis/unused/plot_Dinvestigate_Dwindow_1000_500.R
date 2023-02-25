library(tidyverse)
library(gridExtra)
setwd("~/Desktop/Dstats/Dwindow_1000_500")

#set Z-score variable to identify outlier windows. Higher value is more conservative
ZTHRESH = 4

#generate a data frame which will capture significant D statistic windows
sig_windows <- data.frame()

#Test 1
#first let's look at potential historical introgression between davidsonii and newberryi
############################################################
#read in results
t1 <- read.table("rup_98_new_75_dav_117_localFstats__1000_500.txt", as.is = T, header = T)
t2 <- read.table("rup_98_new_75_dav_87_localFstats__1000_500.txt", as.is = T, header = T)
t3 <- read.table("rup_98_new_80_dav_117_localFstats__1000_500.txt", as.is = T, header = T)
t4 <- read.table("rup_98_new_80_dav_87_localFstats__1000_500.txt", as.is = T, header = T)


#to subset data into different scaffolds ...
#add extra column to each test, specifying which test it is:
t1$testname <- "rup_98_new_75_dav_117"
t2$testname <- "rup_98_new_75_dav_87"
t3$testname <- "rup_98_new_80_dav_117"
t4$testname <- "rup_98_new_80_dav_87"

#generate new data.frame which includes outlier statistics
testlist <- c("t1", "t2", "t3", "t4")

#for loop to generate outlier statistics per scaffold, per sample (test)
for (i in 1:length(testlist)){
  tmptest <- get(testlist[i])
  
  testreplace <- data.frame()
  
  for (j in 1:length(unique(tmptest$chr))){
    tmpscaf <- tmptest %>%
      filter(chr == unique(tmptest$chr)[j]) %>%
      mutate(zscore = (f_dM - mean(f_dM))/sd(f_dM))
    
    testreplace <- rbind(testreplace, tmpscaf)
  }
  
  assign(testlist[i], testreplace)
}


#combine into one large dataset
alltests <- rbind(t1, t2, t3, t4)


#separate by scaffolds
scaflist = unique(alltests$chr)
for (i in 1:length(scaflist)){
  assign(paste("s", i, sep = ""), alltests %>%
           filter(chr == scaflist[i]))
}


#for loop to:
#(1) add outlier windows to a spreadsheet, and 
#(2) generate plots for sliding D statistics across scaffolds

for (i in 1:length(scaflist)){
  #specify in scaffold
  inscaf <- get(paste("s", i, sep = ""))
  
  #identify outliers
  outliers <- rbind(inscaf[which(inscaf$zscore > ZTHRESH),],
                    inscaf[which(inscaf$zscore < -ZTHRESH),])
  sig_windows <- rbind(sig_windows, outliers)
  
  #generate plot
  assign(paste("p", i, sep = ""),
         ggplot(inscaf, aes(x = windowStart/1000000, y = f_dM, color = testname)) +
           facet_wrap(~ testname, ncol = 1) +
           geom_line(linewidth = 0.25) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_point(data=outliers, cex = 1.5, alpha = 0.5) +
           ylab("f_dM") +
           xlab(paste(scaflist[i], " position (Mb)", sep = "")) +
           theme(legend.position = "none"))
}
#write the plot
pdf(paste('davidsonii_x_newberryi_f_dM_window_1000_500_z', ZTHRESH, '.pdf',
          sep = ""))
for (i in 1:length(scaflist)){
  print(get(paste("p", i, sep = "")))
}
dev.off()
############################################################


#Test 2.1
#next let's look at a potential recent introgression between rupicola and newberryi
#here's specifically the rupicola at Mt. Shasta:
############################################################
#read in results
t1 <- read.table("car_91_new_75_rup_86_localFstats__1000_500.txt", as.is = T, header = T)
t2 <- read.table("car_91_new_80_rup_86_localFstats__1000_500.txt", as.is = T, header = T)
t3 <- read.table("car_91_new_85_rup_86_localFstats__1000_500.txt", as.is = T, header = T)


#to subset data into different scaffolds ...
#add extra column to each test, specifying which test it is:
t1$testname <- "car_91_new_75_rup_86"
t2$testname <- "car_91_new_80_rup_86"
t3$testname <- "car_91_new_85_rup_86"

#generate new data.frame which includes outlier statistics
testlist <- c("t1", "t2", "t3")

#for loop to generate outlier statistics per scaffold, per sample (test)
for (i in 1:length(testlist)){
  tmptest <- get(testlist[i])
  
  testreplace <- data.frame()
  
  for (j in 1:length(unique(tmptest$chr))){
    tmpscaf <- tmptest %>%
      filter(chr == unique(tmptest$chr)[j]) %>%
      mutate(zscore = (f_dM - mean(f_dM))/sd(f_dM))
    
    testreplace <- rbind(testreplace, tmpscaf)
  }
  
  assign(testlist[i], testreplace)
}


#combine into one large dataset
alltests <- rbind(t1, t2, t3)


#separate by scaffolds
scaflist = unique(alltests$chr)
for (i in 1:length(scaflist)){
  assign(paste("s", i, sep = ""), alltests %>%
           filter(chr == scaflist[i]))
}


#for loop to:
#(1) add outlier windows to a spreadsheet, and 
#(2) generate plots for sliding D statistics across scaffolds

for (i in 1:length(scaflist)){
  #specify in scaffold
  inscaf <- get(paste("s", i, sep = ""))
  
  #identify outliers
  outliers <- rbind(inscaf[which(inscaf$zscore > ZTHRESH),],
                    inscaf[which(inscaf$zscore < -ZTHRESH),])
  sig_windows <- rbind(sig_windows, outliers)
  
  #generate plot
  assign(paste("p", i, sep = ""),
         ggplot(inscaf, aes(x = windowStart/1000000, y = f_dM, color = testname)) +
           facet_wrap(~ testname, ncol = 1) +
           geom_line(linewidth = 0.25) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_point(data=outliers, cex = 1.5, alpha = 0.5) +
           ylab("f_dM") +
           xlab(paste(scaflist[i], " position (Mb)", sep = "")) +
           theme(legend.position = "none"))
}
#write the plot
pdf(paste('newberryi_x_rupicola_SHASTA_f_dM_window_1000_500_z', ZTHRESH, '.pdf',
          sep = ""))
for (i in 1:length(scaflist)){
  print(get(paste("p", i, sep = "")))
}
dev.off()
############################################################

#Test 2.2
#rupicola and newberryi again, but allopatric pops (not Mt. Shasta)
############################################################
#read in results
t1 <- read.table("car_91_new_80_rup_98_localFstats__1000_500.txt", as.is = T, header = T)
t2 <- read.table("car_91_new_75_rup_98_localFstats__1000_500.txt", as.is = T, header = T)

#to subset data into different scaffolds ...
#add extra column to each test, specifying which test it is:
t1$testname <- "car_91_new_80_rup_98"
t2$testname <- "car_91_new_75_rup_98"

#generate new data.frame which includes outlier statistics
testlist <- c("t1", "t2")

#for loop to generate outlier statistics per scaffold, per sample (test)
for (i in 1:length(testlist)){
  tmptest <- get(testlist[i])
  
  testreplace <- data.frame()
  
  for (j in 1:length(unique(tmptest$chr))){
    tmpscaf <- tmptest %>%
      filter(chr == unique(tmptest$chr)[j]) %>%
      mutate(zscore = (f_dM - mean(f_dM))/sd(f_dM))
    
    testreplace <- rbind(testreplace, tmpscaf)
  }
  
  assign(testlist[i], testreplace)
}


#combine into one large dataset
alltests <- rbind(t1, t2)


#separate by scaffolds
scaflist = unique(alltests$chr)
for (i in 1:length(scaflist)){
  assign(paste("s", i, sep = ""), alltests %>%
           filter(chr == scaflist[i]))
}


#for loop to:
#(1) add outlier windows to a spreadsheet, and 
#(2) generate plots for sliding D statistics across scaffolds

for (i in 1:length(scaflist)){
  #specify in scaffold
  inscaf <- get(paste("s", i, sep = ""))
  
  #identify outliers
  outliers <- rbind(inscaf[which(inscaf$zscore > ZTHRESH),],
                    inscaf[which(inscaf$zscore < -ZTHRESH),])
  sig_windows <- rbind(sig_windows, outliers)
  
  #generate plot
  assign(paste("p", i, sep = ""),
         ggplot(inscaf, aes(x = windowStart/1000000, y = f_dM, color = testname)) +
           facet_wrap(~ testname, ncol = 1) +
           geom_line(linewidth = 0.25) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_point(data=outliers, cex = 1.5, alpha = 0.5) +
           ylab("f_dM") +
           xlab(paste(scaflist[i], " position (Mb)", sep = "")) +
           theme(legend.position = "none"))
}
#write the plot
pdf(paste('newberryi_x_rupicola_NO-SHASTA_f_dM_window_1000_500_z', ZTHRESH, '.pdf',
          sep = ""))
for (i in 1:length(scaflist)){
  print(get(paste("p", i, sep = "")))
}
dev.off()
############################################################


#Test 3
#Crater Lake davidsonii x allopatric rupicola (active hybrid zone)
############################################################
#read in results
t1 <- read.table("dav_118_dav_116_rup_105_localFstats__1000_500.txt", as.is = T, header = T)
t2 <- read.table("dav_118_dav_116_rup_98_localFstats__1000_500.txt", as.is = T, header = T)

#to subset data into different scaffolds ...
#add extra column to each test, specifying which test it is:
t1$testname <- "dav_118_dav_116_rup_105"
t2$testname <- "dav_118_dav_116_rup_98"

#generate new data.frame which includes outlier statistics
testlist <- c("t1", "t2")

#for loop to generate outlier statistics per scaffold, per sample (test)
for (i in 1:length(testlist)){
  tmptest <- get(testlist[i])
  
  testreplace <- data.frame()
  
  for (j in 1:length(unique(tmptest$chr))){
    tmpscaf <- tmptest %>%
      filter(chr == unique(tmptest$chr)[j]) %>%
      mutate(zscore = (f_dM - mean(f_dM))/sd(f_dM))
    
    testreplace <- rbind(testreplace, tmpscaf)
  }
  
  assign(testlist[i], testreplace)
}


#combine into one large dataset
alltests <- rbind(t1, t2)


#separate by scaffolds
scaflist = unique(alltests$chr)
for (i in 1:length(scaflist)){
  assign(paste("s", i, sep = ""), alltests %>%
           filter(chr == scaflist[i]))
}


#for loop to:
#(1) add outlier windows to a spreadsheet, and 
#(2) generate plots for sliding D statistics across scaffolds

for (i in 1:length(scaflist)){
  #specify in scaffold
  inscaf <- get(paste("s", i, sep = ""))
  
  #identify outliers
  outliers <- rbind(inscaf[which(inscaf$zscore > ZTHRESH),],
                    inscaf[which(inscaf$zscore < -ZTHRESH),])
  sig_windows <- rbind(sig_windows, outliers)
  
  #generate plot
  assign(paste("p", i, sep = ""),
         ggplot(inscaf, aes(x = windowStart/1000000, y = f_dM, color = testname)) +
           facet_wrap(~ testname, ncol = 1) +
           geom_line(linewidth = 0.25) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_point(data=outliers, cex = 1.5, alpha = 0.5) +
           ylab("f_dM") +
           xlab(paste(scaflist[i], " position (Mb)", sep = "")) +
           theme(legend.position = "none"))
}
#write the plot
pdf(paste('davidsonii_x_rupicola_CRATERLAKE_f_dM_window_1000_500_z', ZTHRESH, '.pdf',
          sep = ""))
for (i in 1:length(scaflist)){
  print(get(paste("p", i, sep = "")))
}
dev.off()
############################################################


#Test 4 and 5
#fruticosus introgression: fru x rup and fru x car
############################################################
#read in results
t1 <- read.table("dav_118_fru_106_rup_101_localFstats__1000_500.txt", as.is = T, header = T)
t2 <- read.table("dav_118_fru_106_car_28_localFstats__1000_500.txt", as.is = T, header = T)

#to subset data into different scaffolds ...
#add extra column to each test, specifying which test it is:
t1$testname <- "dav_118_fru_106_rup_101"
t2$testname <- "dav_118_fru_106_car_28"

#generate new data.frame which includes outlier statistics
testlist <- c("t1", "t2")

#for loop to generate outlier statistics per scaffold, per sample (test)
for (i in 1:length(testlist)){
  tmptest <- get(testlist[i])
  
  testreplace <- data.frame()
  
  for (j in 1:length(unique(tmptest$chr))){
    tmpscaf <- tmptest %>%
      filter(chr == unique(tmptest$chr)[j]) %>%
      mutate(zscore = (f_dM - mean(f_dM))/sd(f_dM))
    
    testreplace <- rbind(testreplace, tmpscaf)
  }
  
  assign(testlist[i], testreplace)
}


#combine into one large dataset
alltests <- rbind(t1, t2)


#separate by scaffolds
scaflist = unique(alltests$chr)
for (i in 1:length(scaflist)){
  assign(paste("s", i, sep = ""), alltests %>%
           filter(chr == scaflist[i]))
}


#for loop to:
#(1) add outlier windows to a spreadsheet, and 
#(2) generate plots for sliding D statistics across scaffolds

for (i in 1:length(scaflist)){
  #specify in scaffold
  inscaf <- get(paste("s", i, sep = ""))
  
  #identify outliers
  outliers <- rbind(inscaf[which(inscaf$zscore > ZTHRESH),],
                    inscaf[which(inscaf$zscore < -ZTHRESH),])
  sig_windows <- rbind(sig_windows, outliers)
  
  #generate plot
  assign(paste("p", i, sep = ""),
         ggplot(inscaf, aes(x = windowStart/1000000, y = f_dM, color = testname)) +
           facet_wrap(~ testname, ncol = 1) +
           geom_line(linewidth = 0.25) +
           geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
           geom_point(data=outliers, cex = 1.5, alpha = 0.5) +
           ylab("f_dM") +
           xlab(paste(scaflist[i], " position (Mb)", sep = "")) +
           theme(legend.position = "none"))
}
#write the plot
pdf(paste('fruticosus_introgression_f_dM_window_1000_500_z', ZTHRESH, '.pdf',
          sep = ""))
for (i in 1:length(scaflist)){
  print(get(paste("p", i, sep = "")))
}
dev.off()
############################################################


#write the significant windows spreadsheet
write.csv(sig_windows,
          file = paste("significant_windows_1000_500_z", ZTHRESH, ".csv", sep = ""),
          row.names = F)

#write the bedfile which will be used to search for CDS in the windows
bedfile <- sig_windows[, c(1,2,3,8)]
write.table(bedfile,
            file = paste("significant_windows_1000_500_z", ZTHRESH, ".bed", sep = ""),
            row.names = F, col.names = F, quote = F, sep = '\t')

