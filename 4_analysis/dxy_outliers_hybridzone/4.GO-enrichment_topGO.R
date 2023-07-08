#useful:
#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#https://zhiganglu.com/post/topgo-ks-test/

library(tidyverse)
library(topGO)
library(Rgraphviz)

setwd('~/project storage/project_dasanthera_novaseq/results/dxy_outliers_hybridzone_results/')

#genomic background is a tab-delimited file
#first column is gene name
#second column contains a comma-separated list of GO terms
gene2go <- readMappings('topgo_background.tsv', sep = '\t', IDsep = ',')



#next we want to read in our list of genes in outlier windows
#since we have four different groups, the best way to do this is through a for loop
for (j in list.files(path = getwd(), pattern = glob2rx("unique*.bed"))){
  
  #read in gene lists
  outlier_genes <- read.table(j, header = F)
  outlier_genes <- as.character(outlier_genes[,1])
  
  #find where these genes occur in the universe -- created T/F named vector
  genelist <- factor(as.integer(names(gene2go) %in% outlier_genes))
  names(genelist) <- names(gene2go)
  
  #create topGOdata object -- contains list of genes of interest, annotations,
  #and GO hierarchy
  #ontology can be biological processes, molecular function, or cellular component 
  topGO_bp <- new("topGOdata",
                  description = "outlier_fdm_genes",
                  ontology = "BP",
                  nodeSize = 10,
                  allGenes = genelist,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2go)
  
  topGO_mf <- new("topGOdata",
                  description = "outlier_fdm_genes",
                  ontology = "MF",
                  nodeSize = 10,
                  allGenes = genelist,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2go)
  
  topGO_cc <- new("topGOdata",
                  description = "outlier_fdm_genes",
                  ontology = "CC",
                  nodeSize = 10,
                  allGenes = genelist,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2go)
  
  
  #perform fisher's exact tests with a couple different algorithms
  weight01.fisher_bp <- runTest(topGO_bp,
                                algorithm = "weight01",
                                statistic = "fisher")
  
  weight01.fisher_mf <- runTest(topGO_mf,
                                algorithm = "weight01",
                                statistic = "fisher")
  
  weight01.fisher_cc <- runTest(topGO_cc,
                                algorithm = "weight01",
                                statistic = "fisher")
  
  
  #identify significant results at pval
  pval = 0.05
  tbl_bp <- GenTable(topGO_bp, weight01.fisher_bp, topNodes = length(usedGO(object = topGO_bp))) %>%
    filter(result1 < pval)
  tbl_mf <- GenTable(topGO_mf, weight01.fisher_mf, topNodes = length(usedGO(object = topGO_mf))) %>%
    filter(result1 < pval)
  tbl_cc <- GenTable(topGO_cc, weight01.fisher_cc, topNodes = length(usedGO(object = topGO_cc))) %>%
    filter(result1 < pval)
  
  
  #find genes annotated to set of top GO terms
  topgenes_bp <- genesInTerm(topGO_bp, tbl_bp$GO.ID)
  topgenes_mf <- genesInTerm(topGO_mf, tbl_mf$GO.ID)
  topgenes_cc <- genesInTerm(topGO_cc, tbl_cc$GO.ID)
  
  
  #find genes in that set that were part of the initial genes of interest
  #bp
  for(i in 1:length(topgenes_bp)){
    if(i==1){
      genes_in_GO_bp <- vector()
    }
    # find which genes are associated
    myfactor <- which(outlier_genes%in%topgenes_bp[[i]])
    if(length(myfactor) >0){
      genes_in_GO_bp <- append(genes_in_GO_bp, outlier_genes[myfactor])
    }
  }
  
  #mf
  for(i in 1:length(topgenes_mf)){
    if(i==1){
      genes_in_GO_mf <- vector()
    }
    # find which genes are associated
    myfactor <- which(outlier_genes%in%topgenes_mf[[i]])
    if(length(myfactor) >0){
      genes_in_GO_mf <- append(genes_in_GO_mf, outlier_genes[myfactor])
    }
  }
  
  
  #cc
  for(i in 1:length(topgenes_cc)){
    if(i==1){
      genes_in_GO_cc <- vector()
    }
    # find which genes are associated
    myfactor <- which(outlier_genes%in%topgenes_cc[[i]])
    if(length(myfactor) >0){
      genes_in_GO_cc <- append(genes_in_GO_cc, outlier_genes[myfactor])
    }
  }
  
  #create new directory for writing output
  dirbase <- substr(j, nchar(j)-10, nchar(j)-4)
  dirname <- paste("GO_enrichment_output_", dirbase, sep="")
  dir.create(dirname)
  
  
  #write main output files:
  #1. significant GO terms
  write.csv(tbl_bp, file=paste(dirname, "/", dirbase, "_significant_GO_terms_BP.csv", sep = ""))
  write.csv(tbl_mf, file=paste(dirname, "/", dirbase, "_significant_GO_terms_MF.csv", sep = ""))
  write.csv(tbl_cc, file=paste(dirname, "/", dirbase, "_significant_GO_terms_CC.csv", sep = ""))
  
  #2. genes of initial interest with those GO terms included
  write.csv(data.frame(unique(genes_in_GO_bp)), file = paste(dirname, "/", dirbase, "_genes_of_interest_in_significant_GO_terms_BP.csv", sep = ""))
  write.csv(data.frame(unique(genes_in_GO_mf)), file = paste(dirname, "/", dirbase, "_genes_of_interest_in_significant_GO_terms_MF.csv", sep = ""))
  write.csv(data.frame(unique(genes_in_GO_cc)), file = paste(dirname, "/", dirbase, "_genes_of_interest_in_significant_GO_terms_CC.csv", sep = ""))
}












## sidebar
################################################################################
#this is the table of p-values
weight01.fisher_bp@score
min(weight01.fisher_bp@score)

#If I correct p-values for multiple testing with Benjamini Hochberg
#then there are no significant GO enrichment terms
min(p.adjust(weight01.fisher_bp@score, method = "fdr"))
min(p.adjust(weight01.fisher_mf@score, method = "fdr"))
min(p.adjust(weight01.fisher_cc@score, method = "fdr"))



#tbl_bp <- GenTable(topGO_bp, weight01.fisher_bp, topNodes = length(usedGO(object = topGO_bp)))
#tbl_mf <- GenTable(topGO_mf, weight01.fisher_mf, topNodes = length(usedGO(object = topGO_mf)))
#tbl_cc <- GenTable(topGO_cc, weight01.fisher_cc, topNodes = length(usedGO(object = topGO_cc)))


#print pdf of the graph
#printGraph(topGO_bp, weight01.fisher_bp, firstSigNodes = 5, pdfSW = T)
#printGraph(topGO_mf, weight01.fisher_mf, firstSigNodes = 5, pdfSW = T)
#printGraph(topGO_cc, weight01.fisher_cc, firstSigNodes = 5, pdfSW = T)
################################################################################

#find genes annotated with enriched GO terms
weight01.fisher_bp@geneData[1]
weight01.fisher_mf@geneData[1]
weight01.fisher_cc@geneData[1]





