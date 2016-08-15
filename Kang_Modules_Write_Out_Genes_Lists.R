
# Damon Polioudakis
# 2016-08-10
# Write out gene lists from Kang modules as txt files
################################################################################

rm(list = ls())

library(biomaRt)

# Kang
load("../orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
kangExDF = t(datExpr)
rm(datExpr)
# sampleKey is ProcessedKangMetaData.csv
metKangDF <- sampleKey
# metKangDF <- read.csv("../orig.data/InVivoData/ProcessedKangMetaData.csv")
KangAnnotRAW = read.csv("../orig.data/InVivoData/annot.csv", row.names = 1)
KangAnnotnet2 = KangAnnotRAW[which(rownames(KangAnnotRAW) %in% rownames(kangExDF)), ] 
load("../orig.data/InVivoData/net2_power16_cutHeigt0.15.RData")

outDir <- "../analysis/TF_site_enrichment/KangModules_GeneLists"
dir.create(outDir, recursive = TRUE)
################################################################################

# Function to convert list of Entrez IDs to data frame of Ensembl IDs and Gene
# Symbols
ConvertEntrezIdToEnsembl <- function (geneSymList) {
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneDF <- data.frame(geneSymList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id")
    , filters = "entrezgene"
    , values = geneDF
    , mart = ensembl
  )
  ensemblGeneSymDF
}
################################################################################

## Make dataframe of entrez IDs and associated module color
geneColDF <- data.frame(ENTREZ_ID = KangAnnotnet2$ENTREZ_ID
  , MODULE_COLORS = net2$colors)

## Convert entrez IDs to Ensembl IDs and Gene symbols
ensEtzDF <- ConvertEntrezIdToEnsembl(geneColDF$ENTREZ_ID)

## Add ensembl IDs to module color dataframe
geneColDF <- merge(geneColDF, ensEtzDF, by.x = "ENTREZ_ID", by.y = "entrezgene")
geneColDF <- geneColDF[ ,-1]
## Split by modules into list of data frames
geneColLDF <- split(geneColDF, geneColDF$MODULE_COLORS)

## Loop through list, format for VJs pipeline input, and write out as csv
for (module in names(geneColLDF)) {
  df <- data.frame(Genename = geneColLDF[[module]]$hgnc_symbol
    , EnsemblID = geneColLDF[[module]]$ensembl_gene_id)
  df$Sno <- seq(1:nrow(df))
  df <- data.frame(Sno = df$Sno
    , Genename = df$Genename
    , EnsemblID = df$EnsemblID)
  print(head(df))
  write.csv(df, paste0(outDir, "/", module, ".csv")
    , quote = FALSE, row.names = FALSE)
}
################################################################################


