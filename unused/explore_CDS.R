#sort -- 
setwd('~/Desktop/Dwindow_outlier_analysis')
infile <- read.table("combined_unique_genemodels.txt", header = T)

#all I need to do is:
#1. identify all unique upcodes
#2. identify all mRNAs associated with each uniprot code


uplist <- unique(infile[,1])
outmat <- data.frame()
for (i in 1:length(uplist)){
  t <- uplist[i]
  col2 <- infile[which(infile[,1] == t),][,2]
  
  outmat <- rbind(outmat, c(t,  paste(col2, collapse = ",")))
}
colnames(outmat) <- c("uniprot_code", "genemodel")

write.csv(outmat, file = 'non_redundant_genelist.csv')


#this file was copied into the combined output excel file, and then deleted

