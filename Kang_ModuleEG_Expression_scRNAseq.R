
# Damon Polioudakis
# 2016-08-19
# Plot Kang module eigengene expression in Tasic and Zhang scRNA-seq datasets

### TODO
# Convert Tasic mouse to human orthologs

################################################################################

rm(list = ls())

require(WGCNA)
require(ggplot2)
require(biomaRt)
require(reshape2)
require(gridExtra)

options(stringsAsFactors = F)

## Load data and assign variables

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

# Zhang
zhangExDF <- read.table("../../zhang_2016/data/fpkm/mmc3-3.txt", skip = 2)
cNames <- head(read.table("../../zhang_2016/data/fpkm/mmc3-3.txt", sep = "\t", fill = TRUE), 2)
cNames <- cNames[ ,1:63]
cNames <- t(cNames)
cNames[2:5,1] <- cNames[2,1]
cNames[6:9,1] <- cNames[6,1]
cNames[10:15,1] <- cNames[10,1]
cNames[16:27,1] <- cNames[16,1]
cNames[29:33,1] <- cNames[29,1]
cNames[34:36,1] <- cNames[34,1]
cNames[37:38,1] <- cNames[37,1]
cNames[39:42,1] <- cNames[39,1]
cNames[43:44,1] <- cNames[43,1]
cNames[45:48,1] <- cNames[45,1]
cNames[49:50,1] <- cNames[49,1]
cNames[51:56,1] <- cNames[51,1]
cNames[57:58,1] <- cNames[57,1]
cNames[59:60,1] <- cNames[59,1]
cNames[61:63,1] <- cNames[61,1]
# colnames(zhangExDF) <- paste(cNames[ ,1], cNames[ ,2])
colnames(zhangExDF) <- cNames[ ,1]
colnames(zhangExDF)[1] <- "GENE_SYMBOL"

# Tasic
tasicExDF <- read.csv("../../tasic_2016/data/GSE71585_RefSeq_TPM.csv", header = TRUE)
tasicAntDF <- read.csv("../../tasic_2016/metadata/GSE71585_Clustering_Results.csv", header = TRUE)

## Variables
graphCodeTitle <- "Kang_ModuleEG_Expression_scRNAseq.R"
outGraphPfx <- "../analysis/graphs/Kang_ModuleEG_Expression_scRNAseq"

## Output Directory
dir.create("../analysis/graphs", recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

## Function: Convert Ensembl IDs to Gene Symbols
ConvertEntrezToSymbol <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("hgnc_symbol", "entrezgene")
    , filters = "entrezgene"
    , values = moduleGenes
    , mart = ensembl
  )
  ensemblGeneSymDF
}
################################################################################

# Subset to unique entrez Ids
kangEtz <- KangAnnotnet2$ENTREZ_ID[! duplicated(KangAnnotnet2$ENTREZ_ID)]
colors <- net2$colors[! duplicated(KangAnnotnet2$ENTREZ_ID)]

## Zhang MEs
# Convert Zhang Gene Symbols to Entrez and remove duplicates
etzGsymDF <- ConvertEntrezToSymbol(kangEtz)
zhangExDF$GENE_SYMBOL <- toupper(zhangExDF$GENE_SYMBOL)
zhangExDF <- merge(etzGsymDF, zhangExDF, by.x = "hgnc_symbol"
  , by.y = "GENE_SYMBOL")
zhangExDF <- zhangExDF[! duplicated(zhangExDF$entrezgene), ]
# Assign module color to Zhang genes
colorsZhang <- colors[match(kangEtz, as.character(zhangExDF$entrezgene))]
colorsZhang <- colorsZhang[complete.cases(colorsZhang)]
# Calculate MEs
# Human
zhangMEsHs <- moduleEigengenes(t(zhangExDF[ ,c(3:43)]), colorsZhang)$eigengenes
# Mouse
zhangMEsMm <- moduleEigengenes(t(zhangExDF[ ,c(44:64)]), colorsZhang)$eigengenes

## Tasic MEs
# Convert Tasic Gene Symbols to Entrez and remove duplicates
etzGsymDF <- ConvertEntrezToSymbol(kangEtz)
tasicExDF$gene <- toupper(tasicExDF$gene)
tasicExDF <- merge(etzGsymDF, tasicExDF, by.x = "hgnc_symbol"
  , by.y = "gene")
tasicExDF <- tasicExDF[! duplicated(tasicExDF$entrezgene), ]
# Assign module color to Zhang genes
colorsTasic <- colors[match(kangEtz, as.character(tasicExDF$entrezgene))]
colorsTasic <- colorsTasic[complete.cases(colorsTasic)]
# Calculate MEs
# Remove Unclassified cells (were driving variation for some MEs)
keep <- tasicAntDF$broad_type[match(colnames(tasicExDF)[-c(1:2)]
  , gsub("-", "\\.", tasicAntDF$sample_title))] != "Unclassified"
tasicMEs <- moduleEigengenes(t(tasicExDF[-c(1:2)][ ,keep]), colorsTasic)$eigengenes


## Plot

# Zhang Human
df <- data.frame(TYPE = gsub(".\\d+", "", colnames(zhangExDF)[c(3:43)]), zhangMEsHs)
ggDF <- melt(df)
ggDF <- melt(df, id.vars = "TYPE"
  , variable.name = "ME", value.name = "Expression")
ggLDF <- split(ggDF, ggDF$ME)
zhangHsGGL <- lapply(names(ggLDF), function(name) {
  mnGgDF <- aggregate(ggLDF[[name]], list(ggLDF[[name]]$TYPE), mean)
  mnGgDF$TYPE <- mnGgDF$Group.1
  ggplot(ggLDF[[name]], aes(x = TYPE, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.2) +
    geom_point(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "blue") +
    geom_line(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "blue") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Zhang Mouse
df <- data.frame(TYPE = gsub(".\\d+", "", colnames(zhangExDF)[c(44:64)]), zhangMEsMm)
ggDF <- melt(df)
ggDF <- melt(df, id.vars = "TYPE"
  , variable.name = "ME", value.name = "Expression")
ggLDF <- split(ggDF, ggDF$ME)
zhangMmGGL <- lapply(names(ggLDF), function(name) {
  mnGgDF <- aggregate(ggLDF[[name]], list(ggLDF[[name]]$TYPE), mean)
  mnGgDF$TYPE <- mnGgDF$Group.1
  ggplot(ggLDF[[name]], aes(x = TYPE, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.2) +
    geom_point(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "red") +
    geom_line(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "red") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Tasic
# Remove Unclassified cells (were driving variation for some MEs)
type <- tasicAntDF$broad_type[match(colnames(tasicExDF)[-c(1:2)]
  , gsub("-", "\\.", tasicAntDF$sample_title))]
type <- type[type != "Unclassified"]
df <- data.frame(TYPE = type, tasicMEs)
ggDF <- melt(df)
ggDF <- melt(df, id.vars = "TYPE"
  , variable.name = "ME", value.name = "Expression")
ggLDF <- split(ggDF, ggDF$ME)
tasicGGL <- lapply(names(ggLDF), function(name) {
  mnGgDF <- aggregate(ggLDF[[name]], list(ggLDF[[name]]$TYPE), mean)
  mnGgDF$TYPE <- mnGgDF$Group.1
  ggplot(ggLDF[[name]], aes(x = TYPE, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.2) +
    geom_point(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "purple") +
    geom_line(data = mnGgDF, aes(x = TYPE, y = Expression, group = 1), color = "purple") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})


# Combine to set layout of grid.arrange
ggLL <- mapply(list, zhangHsGGL, zhangMmGGL, tasicGGL)

pdf(paste0(outGraphPfx, "Line_Graphs.pdf"), height = 300, width = 10)
do.call("grid.arrange", c(ggLL, ncol = 3))
dev.off()
