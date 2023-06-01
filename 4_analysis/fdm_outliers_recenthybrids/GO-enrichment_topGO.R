library(tidyverse)
library(topGO)
library(Rgraphviz)

setwd('~/project storage/project_dasanthera_novaseq/results/fdm_outlier_analysis/')

#genomic background is a tab-delimited file
#first column is gene name
#second column contains a comma-separated list of GO terms
gene2go <- readMappings('topgo_background.tsv', sep = '\t', IDsep = ',')


##This is where for loop should begin -- loop through each of the window sizes

#next we want to read in our list of genes in outlier windows
outlier_genes <- read.table("unique_mrna_outliers/unique_mRNA_fdm_outliers_1kb.bed",
                            header = F)
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




#list top significant GO terms:
GenTable(topGO_bp, weight01.fisher_bp, topNodes = 10)
GenTable(topGO_mf, weight01.fisher_mf, topNodes = 10)
GenTable(topGO_cc, weight01.fisher_cc, topNodes = 10)




#print pdf of the graph
printGraph(topGOdata, weight01.fisher, firstSigNodes = 5, pdfSW = T)


#find genes annotated with enriched GO terms
weight01.fisher@geneData[1]


#find genes annotated to set of top GO terms
topterms <- sigresults[1:2,1]
topgenes <- genesInTerm(topGOdata, topterms)
topgenes


#correct for multiple testing with Benjamini & Hochberg test 
p.adjust(p, method = "BH", n = length(p))







